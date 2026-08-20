use parry3d::math::{Real, Vector};
use parry3d::shape::{Ball, Capsule, Cone, Cuboid, TriMesh, TriMeshFlags};
use parry3d::transformation::{volume_mesh, VolumeMesh, VolumeMeshParameters};
use std::collections::HashMap;

fn tet_volume(mesh: &VolumeMesh, cell: [u32; 4]) -> Real {
    let [a, b, c, d] = cell.map(|i| mesh.vertices[i as usize]);
    (b - a).cross(c - a).dot(d - a) / 6.0
}

fn total_volume(mesh: &VolumeMesh) -> Real {
    mesh.cells.iter().map(|c| tet_volume(mesh, *c)).sum()
}

/// The faces shared by a single cell: the boundary of the volume mesh.
fn boundary(mesh: &VolumeMesh) -> Vec<[u32; 3]> {
    let mut faces: HashMap<[u32; 3], (usize, [u32; 3])> = HashMap::new();

    for cell in &mesh.cells {
        for face in [
            [cell[0], cell[2], cell[1]],
            [cell[0], cell[1], cell[3]],
            [cell[1], cell[2], cell[3]],
            [cell[0], cell[3], cell[2]],
        ] {
            let mut key = face;
            key.sort_unstable();
            let entry = faces.entry(key).or_insert((0, face));
            entry.0 += 1;
        }
    }

    for (count, _) in faces.values() {
        assert!(*count <= 2, "non-manifold face shared by {count} cells");
    }

    faces
        .values()
        .filter(|(count, _)| *count == 1)
        .map(|(_, face)| *face)
        .collect()
}

/// The extreme dihedral angles of the mesh, in degrees.
fn dihedral_angle_range(mesh: &VolumeMesh) -> (Real, Real) {
    let (mut min, mut max): (Real, Real) = (Real::MAX, 0.0);

    for cell in &mesh.cells {
        let pts = cell.map(|i| mesh.vertices[i as usize]);
        // The outward normal of the face opposite to each vertex.
        let normals: Vec<_> = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|i| pts[i]);
                (b - a).cross(c - a).normalize()
            })
            .collect();

        for i in 0..4 {
            for j in i + 1..4 {
                // The dihedral angle along the edge shared by the two faces.
                let angle = normals[i].dot(normals[j]).clamp(-1.0, 1.0).acos();
                let angle = (core::f64::consts::PI as Real - angle).to_degrees();
                min = min.min(angle);
                max = max.max(angle);
            }
        }
    }

    (min, max)
}

fn check_mesh(mesh: &VolumeMesh, name: &str) {
    assert!(!mesh.cells.is_empty(), "{name}: no cell");

    for cell in &mesh.cells {
        assert!(
            tet_volume(mesh, *cell) > 0.0,
            "{name}: inverted or degenerate cell {cell:?}"
        );
    }

    let (min_angle, max_angle) = dihedral_angle_range(mesh);
    // The angle bounds of the pristine lattice cells (the cover cuts and warps nothing).
    assert!(
        min_angle > 9.05,
        "{name}: minimum dihedral angle {min_angle} is too small"
    );
    assert!(
        max_angle < 160.54,
        "{name}: maximum dihedral angle {max_angle} is too large"
    );
}

#[test]
fn volume_mesh_open_mesh_is_rejected() {
    // A single triangle encloses nothing, so the cover, which needs a sign for its
    // interior, refuses it; the crust is the mode that takes open input.
    let vertices = vec![Vector::ZERO, Vector::X, Vector::Y];
    let indices = vec![[0, 1, 2]];
    assert!(volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(0.1)).is_none());
}

#[test]
fn volume_mesh_rejects_invalid_parameters() {
    let (vertices, indices) = Ball::new(1.0).to_trimesh(10, 10);
    assert!(volume_mesh(&vertices, &indices, &VolumeMeshParameters::new(0.0)).is_none());
    assert!(volume_mesh(&[], &[], &VolumeMeshParameters::new(0.1)).is_none());
}

#[test]
fn volume_mesh_welds_duplicated_vertices() {
    // The same cuboid, with every triangle carrying its own copy of its vertices: closed
    // geometrically, but not topologically.
    let (vertices, indices) = Cuboid::new(Vector::splat(0.5)).to_trimesh();
    let mut split_vertices = Vec::new();
    let split_indices: Vec<[u32; 3]> = indices
        .iter()
        .map(|tri| {
            let first = split_vertices.len() as u32;
            split_vertices.extend(tri.iter().map(|i| vertices[*i as usize]));
            [first, first + 1, first + 2]
        })
        .collect();

    let params = VolumeMeshParameters::new(0.2);
    let mesh = volume_mesh(&split_vertices, &split_indices, &params).unwrap();
    check_mesh(&mesh, "welded cuboid");

    let reference = volume_mesh(&vertices, &indices, &params).unwrap();
    assert_eq!(mesh.cells.len(), reference.cells.len());
}

/// Whether the cell (positively oriented) contains the point, within a relative tolerance.
fn cell_contains(mesh: &VolumeMesh, cell: [u32; 4], pt: Vector) -> bool {
    let [a, b, c, d] = cell.map(|i| mesh.vertices[i as usize]);
    let volume = (b - a).cross(c - a).dot(d - a);
    let tol = -volume.abs() * 1.0e-4;
    (b - pt).cross(c - pt).dot(d - pt) >= tol
        && (pt - a).cross(c - a).dot(d - a) >= tol
        && (b - a).cross(pt - a).dot(d - a) >= tol
        && (b - a).cross(c - a).dot(pt - a) >= tol
}

/// Every vertex, edge midpoint and triangle centroid of the boundary is inside some cell.
fn assert_encloses(mesh: &VolumeMesh, vertices: &[Vector], indices: &[[u32; 3]], name: &str) {
    let samples = vertices
        .iter()
        .copied()
        .chain(indices.iter().flat_map(|idx| {
            let [a, b, c] = idx.map(|i| vertices[i as usize]);
            [
                (a + b) * 0.5,
                (b + c) * 0.5,
                (c + a) * 0.5,
                (a + b + c) / 3.0,
            ]
        }));
    let aabbs: Vec<(Vector, Vector)> = mesh
        .cells
        .iter()
        .map(|cell| {
            let pts = cell.map(|i| mesh.vertices[i as usize]);
            (
                pts.iter().fold(Vector::splat(Real::MAX), |a, b| a.min(*b)),
                pts.iter().fold(Vector::splat(-Real::MAX), |a, b| a.max(*b)),
            )
        })
        .collect();

    'samples: for (k, pt) in samples.enumerate() {
        for (cell, (mins, maxs)) in mesh.cells.iter().zip(&aabbs) {
            if pt.clamp(*mins, *maxs) == pt && cell_contains(mesh, *cell, pt) {
                continue 'samples;
            }
        }
        panic!("{name}: boundary sample {k} at {pt:?} is outside the mesh");
    }
}

/// The cover keeps whole lattice cells, so it contains the shape, subdivided or not, and
/// the unsubdivided cells keep the lattice's own angle bounds (nothing was warped or cut).
#[test]
fn volume_mesh_cover_encloses() {
    use parry3d::transformation::MeshEnclosure;

    let (vertices, indices) = Ball::new(1.0).to_trimesh(20, 20);
    let mut params = VolumeMeshParameters::new(0.25);

    for subdivisions in [0, 1] {
        params.cover_subdivisions = subdivisions;
        params.enclosure = MeshEnclosure::Cover;
        let cover = volume_mesh(&vertices, &indices, &params).unwrap();

        let name = format!("cover ball ({subdivisions} subdivisions)");
        if subdivisions == 0 {
            check_mesh(&cover, &name);
        } else {
            for cell in &cover.cells {
                assert!(tet_volume(&cover, *cell) > 0.0, "inverted cell {cell:?}");
            }
        }
        assert_encloses(&cover, &vertices, &indices, &name);

        // Fatter than the ball it contains, by less than the cell-sized shell around it.
        let expected: Real = 4.0 / 3.0 * core::f64::consts::PI as Real;
        let volume = total_volume(&cover);
        assert!(
            volume > expected && volume < expected * 2.5,
            "cover volume {volume} vs the ball's {expected}"
        );
    }
}

/// The mean distance from the mesh's boundary vertices to the unit sphere.
fn mean_boundary_distance_to_unit_sphere(mesh: &VolumeMesh) -> Real {
    let mut on_boundary = vec![false; mesh.vertices.len()];
    for face in boundary(mesh) {
        for v in face {
            on_boundary[v as usize] = true;
        }
    }
    let (mut sum, mut count) = (0.0, 0);
    for (v, pt) in mesh.vertices.iter().enumerate() {
        if on_boundary[v] {
            sum += (pt.length() - 1.0).abs();
            count += 1;
        }
    }
    sum / count.max(1) as Real
}

/// The shrink-wrap flattens the cover's staircase (its boundary ends up far closer to the
/// surface) without giving up the containment guarantee or inverting a cell.
#[test]
fn volume_mesh_cover_smoothing() {
    use parry3d::transformation::MeshEnclosure;

    let (vertices, indices) = Ball::new(1.0).to_trimesh(20, 20);
    let mut params = VolumeMeshParameters::new(0.25);
    params.enclosure = MeshEnclosure::Cover;
    let raw = volume_mesh(&vertices, &indices, &params).unwrap();
    params.cover_smoothing = 20;
    let smoothed = volume_mesh(&vertices, &indices, &params).unwrap();

    for cell in &smoothed.cells {
        assert!(tet_volume(&smoothed, *cell) > 0.0, "inverted cell {cell:?}");
    }
    assert_encloses(&smoothed, &vertices, &indices, "smoothed cover ball");
    assert!(
        total_volume(&smoothed) < total_volume(&raw),
        "the wrap did not tighten the cover"
    );

    let raw_distance = mean_boundary_distance_to_unit_sphere(&raw);
    let smoothed_distance = mean_boundary_distance_to_unit_sphere(&smoothed);
    println!(
        "cover boundary distance to the sphere: raw {raw_distance:.4}, smoothed {smoothed_distance:.4}"
    );
    assert!(
        smoothed_distance < raw_distance * 0.5,
        "the wrap left the boundary {smoothed_distance} from the sphere against {raw_distance} raw"
    );

    // The demo's configuration: subdivided, smoothed cover.
    let (vertices, indices) = Cone::new(1.5, 1.0).to_trimesh(40);
    let mut params = VolumeMeshParameters::new(0.16);
    params.enclosure = MeshEnclosure::Cover;
    params.cover_subdivisions = 1;
    params.cover_smoothing = 20;
    let cone = volume_mesh(&vertices, &indices, &params).unwrap();
    for cell in &cone.cells {
        assert!(tet_volume(&cone, *cell) > 0.0, "inverted cell {cell:?}");
    }
    assert_encloses(&cone, &vertices, &indices, "smoothed adaptive cover cone");
}

/// Every edge of the mesh's boundary belongs to an even number of boundary faces: a
/// hanging vertex (a crack between subdivision levels) would leave the subdivided side's
/// half-edges unmatched, which a simulation mesh cannot afford.
fn assert_conforming_boundary(mesh: &VolumeMesh, name: &str) {
    let mut edges: HashMap<[u32; 2], u32> = HashMap::new();
    for face in boundary(mesh) {
        for k in 0..3 {
            let (a, b) = (face[k], face[(k + 1) % 3]);
            *edges.entry([a.min(b), a.max(b)]).or_insert(0) += 1;
        }
    }
    for (edge, count) in &edges {
        assert!(
            count % 2 == 0,
            "{name}: boundary edge {edge:?} borders {count} faces"
        );
    }
}

/// Boundary-crossing cells refine below the cell size, outside sub-cells are discarded,
/// and the mesh stays conforming: the staircase shrinks with each subdivision, containment
/// holds throughout, and the shrink-wrap on top hugs tighter than without subdivision.
#[test]
fn volume_mesh_cover_subdivision() {
    use parry3d::transformation::MeshEnclosure;

    let (vertices, indices) = Ball::new(1.0).to_trimesh(16, 16);
    let mut params = VolumeMeshParameters::new(0.35);
    params.enclosure = MeshEnclosure::Cover;

    let mut distances = Vec::new();
    let mut volumes = Vec::new();
    for subdivisions in [0, 1, 2] {
        params.cover_subdivisions = subdivisions;
        let mesh = volume_mesh(&vertices, &indices, &params).unwrap();

        for cell in &mesh.cells {
            assert!(tet_volume(&mesh, *cell) > 0.0, "inverted cell {cell:?}");
        }
        let name = format!("cover ball, {subdivisions} subdivisions");
        assert_encloses(&mesh, &vertices, &indices, &name);
        assert_conforming_boundary(&mesh, &name);

        distances.push(mean_boundary_distance_to_unit_sphere(&mesh));
        volumes.push(total_volume(&mesh));
        println!(
            "{subdivisions} subdivisions: {} cells, boundary distance {:.4}, volume {:.4}",
            mesh.cells.len(),
            distances[subdivisions as usize],
            volumes[subdivisions as usize],
        );
    }
    for s in 1..3 {
        assert!(
            distances[s] < distances[s - 1] * 0.7,
            "subdivision {s} left the staircase at {} against {}",
            distances[s],
            distances[s - 1]
        );
        assert!(volumes[s] < volumes[s - 1]);
    }

    // The wrap on the subdivided cover: its guard scales with the finer boundary cells, so
    // it ends up closer to the surface than the wrap on the unsubdivided one.
    params.cover_smoothing = 20;
    params.cover_subdivisions = 0;
    let smoothed = volume_mesh(&vertices, &indices, &params).unwrap();
    params.cover_subdivisions = 1;
    let subdivided = volume_mesh(&vertices, &indices, &params).unwrap();

    assert_encloses(
        &subdivided,
        &vertices,
        &indices,
        "smoothed subdivided cover",
    );
    assert_conforming_boundary(&subdivided, "smoothed subdivided cover");
    let coarse = mean_boundary_distance_to_unit_sphere(&smoothed);
    let fine = mean_boundary_distance_to_unit_sphere(&subdivided);
    println!("smoothed boundary distance: unsubdivided {coarse:.4}, subdivided {fine:.4}");
    assert!(
        fine < coarse,
        "the subdivided wrap ({fine}) is no tighter than the unsubdivided one ({coarse})"
    );
}

/// The crust covers the surface alone: an open mesh, which every other lattice mode
/// refuses, comes back as a shell of cells hugging its surface, hollow inside, conforming,
/// and the shrink-wrap still applies.
#[test]
fn volume_mesh_crust_covers_open_meshes() {
    use parry3d::transformation::MeshEnclosure;

    // A ball with its cap cut off: boundary edges, so the closed-mesh modes refuse it.
    let (vertices, indices) = Ball::new(1.0).to_trimesh(20, 20);
    let open_indices: Vec<[u32; 3]> = indices
        .iter()
        .copied()
        .filter(|tri| tri.iter().all(|v| vertices[*v as usize].y < 0.8))
        .collect();
    assert!(open_indices.len() < indices.len());

    let mut params = VolumeMeshParameters::new(0.25);
    params.cover_subdivisions = 1;
    assert!(
        volume_mesh(&vertices, &open_indices, &params).is_none(),
        "the cover was expected to refuse the open mesh"
    );

    params.enclosure = MeshEnclosure::Crust;
    params.cover_smoothing = 20;
    let crust = volume_mesh(&vertices, &open_indices, &params).unwrap();

    for cell in &crust.cells {
        assert!(tet_volume(&crust, *cell) > 0.0, "inverted cell {cell:?}");
    }
    assert_conforming_boundary(&crust, "open crust");

    // Every cell hugs the surface; none sits deep inside or far outside.
    for cell in &crust.cells {
        let center = cell
            .map(|i| crust.vertices[i as usize])
            .iter()
            .copied()
            .sum::<Vector>()
            / 4.0;
        let distance = (center.length() - 1.0).abs();
        assert!(
            center.y > 0.7 || distance < 0.5,
            "a crust cell sits {distance} from the surface"
        );
    }

    // The open surface is still enclosed: remap it to its own vertex list first.
    let mut remap = vec![u32::MAX; vertices.len()];
    let mut kept_vertices: Vec<Vector> = Vec::new();
    let kept_indices: Vec<[u32; 3]> = open_indices
        .iter()
        .map(|tri| {
            tri.map(|v| {
                if remap[v as usize] == u32::MAX {
                    remap[v as usize] = kept_vertices.len() as u32;
                    kept_vertices.push(vertices[v as usize]);
                }
                remap[v as usize]
            })
        })
        .collect();
    assert_encloses(&crust, &kept_vertices, &kept_indices, "open crust");
}

/// On a closed mesh the crust is the cover minus its interior: hollow by design.
#[test]
fn volume_mesh_crust_is_hollow() {
    use parry3d::transformation::MeshEnclosure;

    let (vertices, indices) = Ball::new(1.0).to_trimesh(20, 20);
    let mut params = VolumeMeshParameters::new(0.25);
    params.enclosure = MeshEnclosure::Crust;
    let crust = volume_mesh(&vertices, &indices, &params).unwrap();
    params.enclosure = MeshEnclosure::Cover;
    let cover = volume_mesh(&vertices, &indices, &params).unwrap();

    assert!(crust.cells.len() < cover.cells.len());
    for cell in &crust.cells {
        let center = cell
            .map(|i| crust.vertices[i as usize])
            .iter()
            .copied()
            .sum::<Vector>()
            / 4.0;
        assert!(
            center.length() > 0.5,
            "a crust cell sits at depth {}",
            1.0 - center.length()
        );
    }
    assert_encloses(&crust, &vertices, &indices, "closed crust");
}
