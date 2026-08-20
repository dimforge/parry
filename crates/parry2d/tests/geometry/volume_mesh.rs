use parry2d::math::{Real, Vector};
use parry2d::transformation::{volume_mesh, VolumeMesh, VolumeMeshParameters};

/// A counter-clockwise circle, as a closed polyline.
fn disk(radius: Real, num_points: usize) -> (Vec<Vector>, Vec<[u32; 2]>) {
    let vertices = (0..num_points)
        .map(|i| {
            let angle = i as Real / num_points as Real * core::f64::consts::TAU as Real;
            Vector::new(angle.cos(), angle.sin()) * radius
        })
        .collect();
    let indices = (0..num_points as u32)
        .map(|i| [i, (i + 1) % num_points as u32])
        .collect();
    (vertices, indices)
}

fn area(mesh: &VolumeMesh) -> Real {
    mesh.cells
        .iter()
        .map(|cell| {
            let [a, b, c] = cell.map(|i| mesh.vertices[i as usize]);
            (b - a).perp_dot(c - a) / 2.0
        })
        .sum()
}

/// The smallest angle of the mesh, in degrees.
fn min_angle(mesh: &VolumeMesh) -> Real {
    let mut min = Real::MAX;

    for cell in &mesh.cells {
        let pts = cell.map(|i| mesh.vertices[i as usize]);

        for k in 0..3 {
            let (a, b, c) = (pts[k], pts[(k + 1) % 3], pts[(k + 2) % 3]);
            let angle = (b - a)
                .normalize()
                .dot((c - a).normalize())
                .clamp(-1.0, 1.0)
                .acos();
            min = min.min(angle.to_degrees());
        }
    }

    min
}

#[test]
fn volume_mesh_disk() {
    let (vertices, indices) = disk(1.0, 64);
    let mesh = volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(0.15)).unwrap();

    for cell in &mesh.cells {
        let [a, b, c] = cell.map(|i| mesh.vertices[i as usize]);
        assert!(
            (b - a).perp_dot(c - a) > 0.0,
            "inverted or degenerate cell {cell:?}"
        );
    }

    // The boundary is reproduced exactly, so the mesh covers the polygon and nothing else.
    let expected: Real = (0..64)
        .map(|i| {
            let (a, b) = (vertices[i], vertices[(i + 1) % 64]);
            a.perp_dot(b) / 2.0
        })
        .sum();
    assert!(
        (area(&mesh) - expected).abs() < expected * 1.0e-3,
        "disk area {} vs {expected}",
        area(&mesh)
    );

    // The Delaunay refinement targets 30 degrees; it may stop a bit short near the boundary.
    let min = min_angle(&mesh);
    assert!(min > 20.0, "smallest angle {min}");

    // The elements are sized as requested.
    let mean_area = area(&mesh) / mesh.cells.len() as Real;
    assert!(mean_area < 0.15 * 0.15, "mean cell area {mean_area}");
}

#[test]
fn volume_mesh_disk_with_hole() {
    let (mut vertices, mut indices) = disk(1.0, 48);
    let (inner_vertices, inner_indices) = disk(0.4, 24);
    let offset = vertices.len() as u32;
    vertices.extend(inner_vertices);
    indices.extend(inner_indices.iter().map(|e| [e[0] + offset, e[1] + offset]));

    let mesh = volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(0.1)).unwrap();

    // The hole is left empty, whichever way it winds.
    let ring: Real = core::f64::consts::PI as Real * (1.0 - 0.4 * 0.4);
    assert!(
        (area(&mesh) - ring).abs() < ring * 0.02,
        "ring area {} vs {ring}",
        area(&mesh)
    );

    for cell in &mesh.cells {
        let center = cell
            .iter()
            .map(|i| mesh.vertices[*i as usize])
            .sum::<Vector>()
            / 3.0;
        assert!(center.length() > 0.35, "cell inside the hole");
    }
}

/// Asking for smaller elements must not make them worse: the refinement's vertex budget has to
/// follow the requested element size, or it gives up mid-refinement.
#[test]
fn volume_mesh_quality_holds_at_a_fine_resolution() {
    let (vertices, indices) = disk(1.0, 64);

    for cell_size in [0.4, 0.2, 0.1, 0.05] {
        let mesh = volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(cell_size)).unwrap();
        let min = min_angle(&mesh);
        assert!(min > 20.0, "smallest angle {min} at cell size {cell_size}");

        let mean_area = area(&mesh) / mesh.cells.len() as Real;
        assert!(
            mean_area < cell_size * cell_size,
            "mean cell area {mean_area} at cell size {cell_size}"
        );
    }
}

#[test]
fn volume_mesh_rejects_invalid_parameters() {
    let (vertices, indices) = disk(1.0, 16);
    assert!(volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(0.0)).is_none());
    assert!(volume_mesh(&[], &[], &VolumeMeshParameters::new(0.1)).is_none());
}
