//! Visualises the polyline orientation bug and its fix, side by side.
//!
//! A [`Polyline`] is double-sided by default: a contact normal can point to whichever side of a
//! segment the other shape happens to be nearest, so a body that crosses the line is pushed deeper
//! through it instead of back out.  An *oriented* polyline (built with
//! [`PolylineFlags::ORIENTED`]) clamps every contact normal into the outward cone, so a
//! one-sided boundary pushes a body back to the play side no matter which side it is on.
//!
//! The same boundaries (a convex peak with solid below, and a box with solid inside whose 90-degree
//! corners are where ghost collisions are worst) and the same circle path are drawn twice: the upper
//! panel uses double-sided polylines (the bug), the lower panel uses oriented ones (the fix).  The
//! circle's motion has three phases: a free lissajous roam, a guided trace that hugs each segment and
//! presses into every corner, then a guided trace that only traces the outside (never entering).
//! Each frame computes the contact manifold
//! (the same path the physics engine uses; this is where pseudo-normals apply, not the simpler
//! `query::contact`) and draws the actual contact normal and where the solver would push the circle.
//!
//! The test circle turns magenta when its center is *inside* a boundary, with a line to the nearest
//! surface of each. Only the oriented (lower) polylines have an interior, so only they light up.
//!
//! ```bash
//! cargo run -p parry2d --example oriented_polyline2d
//! ```

mod utils2d;

use core::f32::consts::{PI, TAU};

use kiss3d::prelude::*;
use parry2d::math::{Pose, Rotation};
use parry2d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
    PointQuery,
};
use parry2d::shape::{Ball, Polyline, PolylineFlags, Shape};
use utils2d::{draw_circle, draw_line_2d, lissajous_2d};

/// The test circle's fill when its center is inside an oriented boundary.
const INSIDE_COLOR: Color = Color::new(1.0, 0.30, 0.80, 1.0);

/// One stop on the guided trace: where to be, how long to take getting there, and how long to dwell.
struct Waypoint {
    pos: Vec2,
    travel: f32,
    dwell: f32,
}

#[kiss3d::main]
async fn main() {
    print_legend();

    let mut window =
        Window::new("oriented_polyline2d -- double-sided (top) vs oriented (bottom)").await;
    let mut camera = PanZoomCamera2d::new(Vec2::new(0.0, 2.0), 2.0);
    let mut scene = SceneNode2d::empty();
    let font = Font::default();

    // Convex peak, wound right-to-left so outward (right-hand) normals point up: solid below, play above.
    let peak = Polyline::new(
        vec![
            Vec2::new(96.0, -36.0),
            Vec2::new(54.0, 40.0),
            Vec2::new(12.0, -36.0),
        ],
        Some(vec![[0, 1], [1, 2]]),
    );

    // Box, wound counter-clockwise so the solid is inside; its square corners are the issue's case.
    let boxed = Polyline::new(
        vec![
            Vec2::new(-12.0, -36.0),
            Vec2::new(-12.0, 36.0),
            Vec2::new(-96.0, 36.0),
            Vec2::new(-96.0, -36.0),
        ],
        Some(vec![[0, 1], [1, 2], [2, 3], [3, 0]]),
    );

    // The two panels share geometry; only the lower one is told its outward side.
    let plain = [peak, boxed];
    let mut oriented = plain.clone();
    for boundary in &mut oriented {
        boundary.set_flags(PolylineFlags::ORIENTED);
    }

    let ball = Ball::new(16.0);

    // Three phases: free roam, a corner-pressing trace, then an outside-only trace.
    let trace_in = build_trace(&plain, ball.radius, -0.5 * ball.radius);
    let trace_out = build_trace(&plain, ball.radius, ball.radius);
    let trace_in_duration: f32 = trace_in
        .iter()
        .map(|waypoint| waypoint.travel + waypoint.dwell)
        .sum();
    let trace_out_duration: f32 = trace_out
        .iter()
        .map(|waypoint| waypoint.travel + waypoint.dwell)
        .sum();
    let lissajous_duration = 8.0;
    let cycle = lissajous_duration + trace_in_duration + trace_out_duration;

    let dispatcher = DefaultQueryDispatcher;
    let panel = Vec2::new(0.0, 81.0); // upper panel sits at +panel, lower at -panel

    let start_time = web_time::Instant::now();

    while window.render_2d(&mut scene, &mut camera).await {
        let elapsed = start_time.elapsed().as_secs_f32();

        // The same path feeds both panels.
        let phase = if cycle > 0.0 { elapsed % cycle } else { 0.0 };
        let local_pos = if phase < lissajous_duration {
            lissajous_2d(phase * 0.35) * 120.0
        } else if phase < lissajous_duration + trace_in_duration {
            trace_position(&trace_in, phase - lissajous_duration)
        } else {
            trace_position(&trace_out, phase - lissajous_duration - trace_in_duration)
        };

        draw_panel(&mut window, &dispatcher, &plain, &ball, local_pos, panel);
        draw_panel(
            &mut window,
            &dispatcher,
            &oriented,
            &ball,
            local_pos,
            -panel,
        );

        window.draw_text(
            "upper: double-sided polyline (unoriented)",
            Vec2::new(8.0, 6.0),
            34.0,
            &font,
            WHITE,
        );
        window.draw_text(
            "lower: oriented polyline (with pseudo-normals)",
            Vec2::new(8.0, 42.0),
            34.0,
            &font,
            WHITE,
        );
        window.draw_text(
            "magenta line + circle: center is inside that boundary (oriented panel only)",
            Vec2::new(8.0, 78.0),
            34.0,
            &font,
            INSIDE_COLOR,
        );
    }
}

/// Draws one panel (boundaries, circle, contact resolution), translated by `offset`.
fn draw_panel(
    window: &mut Window,
    dispatcher: &DefaultQueryDispatcher,
    boundaries: &[Polyline],
    ball: &Ball,
    local_pos: Vec2,
    offset: Vec2,
) {
    let radius = ball.radius;

    for boundary in boundaries {
        draw_boundary(window, boundary, offset);
        draw_cones(window, boundary, offset);
    }

    // Draw a line to the nearest surface of each boundary the center is inside (a point can be inside
    // several), then tint the circle. Only oriented polylines report an interior.
    let mut inside_any = false;
    for boundary in boundaries {
        let projection = boundary.project_local_point(local_pos, false);
        if projection.is_inside {
            inside_any = true;
            draw_line_2d(
                window,
                local_pos + offset,
                projection.point + offset,
                INSIDE_COLOR,
            );
        }
    }
    let circle_color = if inside_any { INSIDE_COLOR } else { YELLOW };
    draw_circle(window, local_pos + offset, radius, circle_color);

    for boundary in boundaries {
        draw_resolution(
            window,
            dispatcher,
            boundary,
            ball,
            local_pos,
            offset,
            |w, p| {
                draw_circle(w, p, radius, LIME);
            },
        );
    }
}

/// Whether the polyline is a closed loop (last segment's end meets the first's start).
fn is_closed(polyline: &Polyline) -> bool {
    let num = polyline.num_segments();
    num > 0 && (polyline.segment((num - 1) as u32).b - polyline.segment(0).a).length() < 1.0e-3
}

/// Builds a trace hugging each segment; at each corner it moves `corner_offset` along the outward
/// bisector (negative dips inside, positive stays outside) and dwells.
fn build_trace(boundaries: &[Polyline], radius: f32, corner_offset: f32) -> Vec<Waypoint> {
    const CORNER_DWELL: f32 = 0.8;

    let mut waypoints: Vec<Waypoint> = Vec::new();

    for boundary in boundaries {
        let num = boundary.num_segments();
        if num == 0 {
            continue;
        }

        // Closed loops (the box) wrap last->first; open chains (the peak) don't, so their ends aren't corners.
        let closed = is_closed(boundary);

        let normals: Vec<Vec2> = (0..num as u32)
            .map(|i| {
                let segment = boundary.segment(i);
                let direction = (segment.b - segment.a).normalize_or_zero();
                Vec2::new(direction.y, -direction.x)
            })
            .collect();

        for i in 0..num {
            let segment = boundary.segment(i as u32);
            let normal = normals[i];

            push_waypoint(&mut waypoints, segment.a + normal * radius, false, 0.0);
            push_waypoint(&mut waypoints, segment.b + normal * radius, false, 0.0);

            let next = if i + 1 < num {
                Some(i + 1)
            } else if closed {
                Some(0)
            } else {
                None
            };

            if let Some(j) = next {
                let bisector = (normal + normals[j]).normalize_or_zero();
                push_waypoint(
                    &mut waypoints,
                    segment.b + bisector * corner_offset,
                    true,
                    CORNER_DWELL,
                );
            }
        }
    }

    waypoints
}

/// Appends a waypoint, timing the move from the previous: a steady glide, or a slow crawl if `slow`.
fn push_waypoint(waypoints: &mut Vec<Waypoint>, pos: Vec2, slow: bool, dwell: f32) {
    const GLIDE_SPEED: f32 = 60.0; // pixels per second
    const APPROACH_TRAVEL: f32 = 1.4; // seconds to cross into a corner

    let travel = match waypoints.last() {
        None => 0.0,
        Some(_) if slow => APPROACH_TRAVEL,
        Some(previous) => ((pos - previous.pos).length() / GLIDE_SPEED).max(0.05),
    };

    waypoints.push(Waypoint { pos, travel, dwell });
}

/// Position along the guided trace at local time `t` (seconds since the trace phase started).
fn trace_position(waypoints: &[Waypoint], mut t: f32) -> Vec2 {
    let Some(first) = waypoints.first() else {
        return Vec2::ZERO;
    };

    let mut current = first.pos;
    for waypoint in waypoints {
        if t < waypoint.travel {
            let fraction = if waypoint.travel > 0.0 {
                t / waypoint.travel
            } else {
                1.0
            };
            return current.lerp(waypoint.pos, fraction);
        }
        t -= waypoint.travel;

        if t < waypoint.dwell {
            return waypoint.pos;
        }
        t -= waypoint.dwell;

        current = waypoint.pos;
    }

    current
}

/// Draws a boundary (white) with a gray outward-reference arrow per segment, translated by `offset`.
fn draw_boundary(window: &mut Window, polyline: &Polyline, offset: Vec2) {
    for i in 0..polyline.num_segments() as u32 {
        let segment = polyline.segment(i);
        let a = segment.a + offset;
        let b = segment.b + offset;
        draw_line_2d(window, a, b, WHITE);

        let direction = (segment.b - segment.a).normalize_or_zero();
        // Right-hand normal: the outward/play side for this winding.
        let outward = Vec2::new(direction.y, -direction.x);
        draw_arrow(window, (a + b) * 0.5, outward, 22.0, GRAY);
    }
}

/// Draws the valid-normal cone at each corner (the wedge between adjacent face normals) in the
/// oriented panel; a double-sided polyline returns `None` here, so the upper panel gets no cone.
fn draw_cones(window: &mut Window, polyline: &Polyline, offset: Vec2) {
    const CONE_COLOR: Color = Color::new(0.30, 0.60, 1.0, 1.0);

    let num = polyline.num_segments();
    if num == 0 {
        return;
    }

    let closed = is_closed(polyline);

    for i in 0..num {
        let next = if i + 1 < num {
            i + 1
        } else if closed {
            0
        } else {
            continue;
        };

        // The vertex pseudo-normal bisects the two adjacent faces, so the clamp's cone is that wedge.
        let (Some(here), Some(beyond)) = (
            polyline.segment_normal_constraints(i as u32),
            polyline.segment_normal_constraints(next as u32),
        ) else {
            return;
        };

        let vertex = polyline.segment(i as u32).b + offset;
        draw_cone(window, vertex, here.face, beyond.face, CONE_COLOR);
    }
}

/// Draws a fan filling the angular wedge from `dir_a` to `dir_b` (the shorter way) at `vertex`.
fn draw_cone(window: &mut Window, vertex: Vec2, dir_a: Vec2, dir_b: Vec2, color: Color) {
    const LENGTH: f32 = 30.0;
    const RAYS: usize = 6;

    let angle_a = dir_a.y.atan2(dir_a.x);
    let mut delta = dir_b.y.atan2(dir_b.x) - angle_a;
    while delta > PI {
        delta -= TAU;
    }
    while delta < -PI {
        delta += TAU;
    }

    let mut previous: Option<Vec2> = None;
    for k in 0..=RAYS {
        let angle = angle_a + delta * (k as f32 / RAYS as f32);
        let tip = vertex + Vec2::new(angle.cos(), angle.sin()) * LENGTH;
        draw_line_2d(window, vertex, tip, color);
        if let Some(prev) = previous {
            draw_line_2d(window, prev, tip, color);
        }
        previous = Some(tip);
    }
}

/// Draws the deepest contact's normal (red) and, via `draw_ghost`, the resolved position. The manifold
/// is computed in the boundary's local frame; only the drawing is translated by `offset`.
fn draw_resolution(
    window: &mut Window,
    dispatcher: &DefaultQueryDispatcher,
    polyline: &Polyline,
    shape: &dyn Shape,
    local_pos: Vec2,
    offset: Vec2,
    draw_ghost: impl Fn(&mut Window, Vec2),
) {
    const PREDICTION: f32 = 8.0;

    let pos12 = Pose::from_parts(local_pos, Rotation::identity());
    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace: Option<ContactManifoldsWorkspace> = None;
    // Polyline vs ball is always supported, so the unsupported-pair error cannot happen here.
    let _ = dispatcher.contact_manifolds(
        &pos12,
        polyline,
        shape,
        PREDICTION,
        &mut manifolds,
        &mut workspace,
    );

    let mut deepest: Option<(Vec2, Vec2, f32)> = None;
    for manifold in &manifolds {
        let normal = manifold.local_n1;
        for contact in &manifold.points {
            if deepest.is_none_or(|(_, _, best)| contact.dist < best) {
                deepest = Some((contact.local_p1, normal, contact.dist));
            }
        }
    }

    if let Some((point, normal, distance)) = deepest {
        draw_arrow(window, point + offset, normal, 45.0, RED);

        let penetration = (-distance).max(0.0);
        if penetration > 0.0 {
            draw_ghost(window, local_pos + normal * penetration + offset);
        }
    }
}

/// Draws an arrow from `base` along `direction`, scaled to `length`.
fn draw_arrow(window: &mut Window, base: Vec2, direction: Vec2, length: f32, color: Color) {
    let Some(direction) = direction.try_normalize() else {
        return;
    };

    let tip = base + direction * length;
    draw_line_2d(window, base, tip, color);

    let back = -direction * (length * 0.28);
    let side = Vec2::new(-direction.y, direction.x) * (length * 0.16);
    draw_line_2d(window, tip, tip + back + side, color);
    draw_line_2d(window, tip, tip + back - side, color);
}

fn print_legend() {
    println!("oriented_polyline2d orientation demo (side-by-side)");
    println!("  UPPER panel: double-sided polyline -- the bug");
    println!("  LOWER panel: oriented polyline -- the fix");
    println!("both panels run the same circle path; scroll to zoom, drag to pan.");
    println!("legend:");
    println!("  white  boundary polylines: a peak (solid below) and a box (solid inside)");
    println!("  gray   reference outward (play-side) normal per segment");
    println!("  blue   valid-normal cone the oriented polyline allows at each corner (lower only)");
    println!("  yellow test circle (center outside every boundary)");
    println!("  magenta circle + line to each oriented boundary whose interior holds the center (lower only)");
    println!("  red    actual contact normal -- the direction the solver pushes the circle");
    println!("  green  where the circle is pushed out to");
    println!("phases: free roam, then a corner-pressing trace, then an outside-only trace.");
    println!(
        "watch the corners: the upper panel has no cone, so red/green flip to the wrong side; in"
    );
    println!(
        "the lower panel the red normal stays inside the blue cone, pushing the circle back out."
    );
}
