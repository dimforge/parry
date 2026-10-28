use parry3d::math::{Pose, Rotation, Vector};
use parry3d::query::details::contact_manifold_pfm_pfm_shapes;
use parry3d::query::ContactManifold;
use parry3d::shape::{Cuboid, Cylinder};

// A wheel (cylinder, radius 0.1, half-height 0.05) spinning at ~43 rad/s and rolling on a
// 200x200x0.5 ground cuboid. The poses are consecutive simulation steps, so the manifold's
// cached normal is carried from one step to the next, as in a physics engine.
//
// The wheel hovers ~1e-4 above the ground. With the ground's far corners, the CSO support
// points reach ~144, and GJK stalled with a lower bound of -3.3e-2 on the penetration. That
// was accepted as "touching" because the tolerance scaled with `sqrt(eps) * support_scale`
// (~0.15), more than the whole wheel, and the stale `old_dir` was returned as the contact
// normal: (-0.51, 0.69, 0.51) instead of (0, 1, 0).
#[test]
fn cylinder_rolling_on_large_cuboid_keeps_a_vertical_normal() {
    // (body translation, body rotation as xyzw) for 9 consecutive steps.
    const POSES: [([f32; 3], [f32; 4]); 9] = [
        (
            [2.1498301029205322, 0.12448107451200485, 0.3996926546096802],
            [
                0.2489599585533142,
                -0.6622065305709839,
                0.24911077320575714,
                0.6613963842391968,
            ],
        ),
        (
            [2.1495745182037354, 0.12121988832950592, 0.3993242681026459],
            [
                0.46536311507225037,
                -0.5327503085136414,
                0.46630358695983887,
                0.5312016010284424,
            ],
        ),
        (
            [2.149324655532837, 0.11523353308439255, 0.39895594120025635],
            [
                0.6219001412391663,
                -0.335325688123703,
                0.6241053342819214,
                0.333600789308548,
            ],
        ),
        (
            [2.149080514907837, 0.10652196407318115, 0.398587703704834],
            [
                0.698938250541687,
                -0.0951799601316452,
                0.7025531530380249,
                0.09405004233121872,
            ],
        ),
        (
            [2.1422312259674072, 0.10064999759197235, 0.3972579836845398],
            [
                0.6860884428024292,
                0.16269077360630035,
                0.6905015707015991,
                -0.16131402552127838,
            ],
        ),
        (
            [2.1318469047546387, 0.10105381906032562, 0.3973941504955292],
            [
                0.5856216549873352,
                0.39615002274513245,
                0.5881026387214661,
                -0.39274439215660095,
            ],
        ),
        (
            [2.120213747024536, 0.10016228258609772, 0.3976336121559143],
            [
                0.409687876701355,
                0.5779147148132324,
                0.4103825092315674,
                -0.5742447972297668,
            ],
        ),
        (
            [2.105894088745117, 0.1001090407371521, 0.3976670503616333],
            [
                0.1809343695640564,
                0.6852918267250061,
                0.1802733689546585,
                -0.6820118427276611,
            ],
        ),
        (
            [2.0811259746551514, 0.11391204595565796, 0.40611711144447327],
            [
                -0.03625503554940224,
                0.7045886516571045,
                -0.047442466020584106,
                -0.7070997357368469,
            ],
        ),
    ];
    const PREDICTION: f32 = 0.02;

    let ground = Cuboid::new(Vector::new(100.0, 0.25, 100.0));
    let ground_pose = Pose::from_translation(Vector::new(0.0, -0.25, 0.0));
    let wheel = Cylinder::new(0.05, 0.1);
    // The cylinder's axis (local Y) is laid along the body's local X.
    let axle_to_body_x = Rotation::from_xyzw(
        0.0,
        0.0,
        core::f32::consts::FRAC_1_SQRT_2,
        core::f32::consts::FRAC_1_SQRT_2,
    );

    let mut manifold = ContactManifold::<(), ()>::new();
    for (step, (t, q)) in POSES.iter().enumerate() {
        let body = Pose::from_parts(
            Vector::new(t[0], t[1], t[2]),
            Rotation::from_xyzw(q[0], q[1], q[2], q[3]),
        );
        let wheel_pose = body * Pose::from_parts(Vector::ZERO, axle_to_body_x);
        let pos12 = ground_pose.inv_mul(&wheel_pose);

        contact_manifold_pfm_pfm_shapes(
            &pos12,
            &ground,
            &wheel,
            None,
            None,
            PREDICTION,
            &mut manifold,
        );

        assert!(!manifold.points.is_empty(), "step {step}: no contact");
        assert!(
            manifold.local_n1.y > 0.99,
            "step {step}: contact normal {:?} is not vertical",
            manifold.local_n1
        );
    }
}
