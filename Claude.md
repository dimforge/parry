# Parry: 2D and 3D Collision Detection Library

## Overview

**Parry** is a comprehensive 2D and 3D geometric and collision detection library written in Rust by [Dimforge](https://dimforge.com). It provides high-performance geometric queries, collision detection, and spatial partitioning for game engines, physics simulations, robotics, and other applications requiring geometric computations.

**Documentation Status**: All public functions now have comprehensive, beginner-friendly documentation with 300+ working doc-test examples.

### Key Facts

- **Language**: Rust
- **Organization**: Dimforge
- **License**: Apache-2.0
- **Current Version**: 0.25.0
- **Lines of Code**: ~48,000 lines
- **Documentation**: ~15,000 lines of enhanced documentation
- **Doc-Test Examples**: 300+ tested examples
- **Repository**: https://github.com/dimforge/parry
- **Official Documentation**:
  - 2D: http://docs.rs/parry2d
  - 3D: http://docs.rs/parry3d
  - User Guide: http://parry.rs

## Project Architecture

### Unique Structure

Parry uses an unusual architecture to **share the same codebase** between 2D and 3D versions:

1. **Single Source Directory** (`src/`): Contains all `.rs` source code
2. **Multiple Crate Directories** (`crates/`): Contains separate `Cargo.toml` files for each variant
3. **Conditional Compilation**: Uses cargo features (`dim2`, `dim3`, `f32`, `f64`) to compile the appropriate code

This clever design allows:
- Single implementation for most algorithms
- Type-safe dimension handling via generics
- Minimal code duplication
- Consistent API across dimensions

### Available Crates

The library is distributed as **four separate crates**:

- `parry2d` - 2D collision detection with f32 precision
- `parry3d` - 3D collision detection with f32 precision
- `parry2d-f64` - 2D collision detection with f64 precision
- `parry3d-f64` - 3D collision detection with f64 precision

## Core Modules

### 1. `shape` - Geometric Shapes (25+ types)

Defines all geometric primitives and composite shapes supported by Parry.

#### Basic Shapes (Both 2D and 3D)

**`Ball`** - Sphere (3D) or circle (2D)
- Simplest shape, defined by radius only
- Extremely efficient collision detection
- Perfect for projectiles, particles, bounding volumes
- Can scale to ellipse/ellipsoid (approximated as mesh)

**`Cuboid`** - Box (3D) or rectangle (2D)
- Defined by half-extents on each axis
- Axis-aligned in local space
- Supports non-uniform scaling
- Ideal for walls, platforms, containers

**`Capsule`** - Pill shape (cylinder with hemispherical caps)
- Defined by central segment + radius
- Completely smooth surface (no edges)
- Perfect for characters, elongated objects
- Better than cylinders for dynamic objects

**`Segment`** - Line segment between two points
- 1-dimensional (no thickness)
- Building block for polylines and meshes
- Zero volume, zero mass
- Use Capsule if thickness is needed

**`Triangle`** - Three-vertex polygon
- Fundamental mesh building block
- Has area (2D) but no volume (3D)
- Counter-clockwise winding convention
- Supports barycentric coordinates

**`HalfSpace`** - Infinite half-plane
- Divides space into inside/outside regions
- Defined by a normal vector
- Perfect for ground planes, walls, boundaries
- Very efficient (infinite extent)

#### 3D-Only Shapes

**`Cone`** - Tapered shape with circular base
- Apex points upward (+Y axis)
- Flat circular base
- Use cases: traffic cones, funnels, projectile tips

**`Cylinder`** - Circular cross-section, flat caps
- Axis aligned with Y axis
- Flat top and bottom (sharp edges at rims)
- Use cases: pillars, cans, wheels, pipes

**`ConvexPolyhedron`** - Arbitrary 3D convex mesh
- Created from vertices and triangle indices
- Rich topology (vertices, edges, faces, adjacency)
- Coplanar triangles merged into polygonal faces
- Supports complex convex shapes

**`Tetrahedron`** - Four-vertex polyhedron
- Simplest 3D polyhedron
- Has volume and inertia
- Building block for 3D meshes
- Useful for volumetric computations

#### 2D-Only Shapes

**`ConvexPolygon`** - Arbitrary 2D convex polygon
- Created from vertices or via convex hull
- Counter-clockwise winding required
- Stores vertices and edge normals
- Can be offset (dilated/eroded)

#### Composite Shapes (Require `alloc` feature)

**`Compound`** - Collection of sub-shapes with positions
- Combines multiple shapes into one
- Internal BVH for efficient queries
- Perfect for complex objects (tables, robots, vehicles)
- Can decompose triangle meshes into convex pieces

**`TriMesh`** - Triangle mesh (2D or 3D)
- Large triangle collections
- Internal BVH acceleration
- Optional topology and pseudo-normals
- Flags for internal edge handling
- Supports mesh-mesh boolean operations (with `spade` feature)

**`Polyline`** - Connected line segments
- 1-dimensional path representation
- Automatic or custom connectivity
- Internal BVH for ray casting
- Use cases: paths, boundaries, wireframes

**`HeightField`** - Terrain height map
- 2D: Array of heights (terrain profile)
- 3D: Grid of heights (terrain surface)
- Very efficient for large terrains
- Supports cell removal (holes, water)

**`Voxels`** - Voxel grid representation
- Sparse chunk-based storage (as of v0.25.0)
- Automatic internal edge resolution
- Dynamic modification support
- Perfect for destructible terrain, Minecraft-like worlds
- Memory efficient for large grids

#### Shape Utilities

**`SharedShape`** - Arc-wrapped shape for sharing
- Enables multiple colliders sharing the same geometry
- Copy-on-write via `make_mut()`
- Cheap cloning (reference counting)
- Type erasure for heterogeneous collections

**`RoundShape<S>`** - Adds border radius to any shape
- Creates "padded" version of a shape
- Implements Minkowski sum with a ball
- Softer collisions, better numerical stability
- Examples: `RoundCuboid`, `RoundTriangle`, `RoundCylinder`

**Supporting Types**:
- `FeatureId` - Identifies vertices, edges, faces on shapes
- `PackedFeatureId` - Memory-efficient feature ID (bit-packed)
- `PolygonalFeature` - Polygonal feature with vertices and normals

#### Key Traits

**`Shape`** - Main trait all shapes implement
- AABB and bounding sphere computation
- Mass properties calculation
- Ray casting and point queries
- Shape cloning, scaling, and type identification
- Feature normal queries

**`SupportMap`** - Enables GJK/EPA algorithms
- Provides support point in any direction
- Core interface for convex collision detection
- Extremely efficient (constant time for most shapes)
- Used by: Ball, Cuboid, Capsule, Segment, Triangle, Cylinder, Cone, ConvexPolygon, ConvexPolyhedron

**`CompositeShape`** - For shapes with multiple parts
- Enables BVH-accelerated queries
- Implemented by: Compound, TriMesh, Polyline, HeightField

**`PolygonalFeatureMap`** - For shapes with polygonal features
- Provides vertices, edges, and normals
- Used for accurate contact manifold generation

### 2. `query` - Geometric Queries (50+ functions)

Non-persistent geometric queries between shapes.

#### Main Query Functions

**`distance(pos1, shape1, pos2, shape2)`**
- Returns minimum distance between shapes
- 0.0 if touching or penetrating
- Very fast for simple shapes (Ball-Ball)
- Uses GJK for convex shapes

**`closest_points(pos1, shape1, pos2, shape2, max_dist)`**
- Returns closest point pair
- Three outcomes: `Intersecting`, `WithinMargin`, `Disjoint`
- Useful for proximity detection
- More information than `distance()` alone

**`contact(pos1, shape1, pos2, shape2, prediction)`**
- Returns contact point, normal, and penetration depth
- Supports prediction distance for CCD
- Essential for physics simulation
- Uses GJK+EPA for convex shapes

**`intersection_test(pos1, shape1, pos2, shape2)`**
- Fastest query - boolean only
- Early-exit optimization
- Perfect for broad-phase filtering
- Use for trigger volumes, overlap detection

**`cast_shapes(pos1, vel1, shape1, pos2, vel2, shape2, options)`**
- Linear shape casting (swept collision)
- Finds time of first impact during linear motion
- Essential for continuous collision detection (CCD)
- Prevents tunneling at high speeds
- Returns: time, witness points, normals

**`cast_shapes_nonlinear(motion1, shape1, motion2, shape2, ...)`**
- Shape casting with rotation (6-DOF motion)
- Handles spinning, tumbling objects
- More expensive than linear version
- Necessary for rotating objects at high speeds

#### Query Traits

**`RayCast`** - Ray casting against shapes
- `cast_ray()` - Returns time of impact and normal
- `cast_ray_and_get_normal()` - Full intersection data
- `intersects_ray()` - Boolean check
- Solid parameter: treat interior differently

**`PointQuery`** - Point projection onto shapes
- `project_point()` - Find closest surface point
- `distance_to_point()` - Distance calculation
- `contains_point()` - Interior test
- `project_point_and_get_feature()` - With feature ID

**Supporting Types**:
- `Ray` - Origin + direction for ray casting
- `RayIntersection` - Result with time, normal, feature
- `PointProjection` - Result with point and is_inside flag
- `ClosestPoints` - Enum for closest point results

#### Collision Detection Algorithms

**GJK (Gilbert-Johnson-Keerthi)**
- Computes distance between convex shapes
- Iterative simplex refinement
- Works on any `SupportMap` shapes
- Very fast convergence (usually <10 iterations)
- Cannot compute penetration depth (use EPA)

**EPA (Expanding Polytope Algorithm)**
- Computes penetration depth for overlapping convex shapes
- Continuation of GJK when shapes penetrate
- Iteratively expands polytope to find minimum translation
- Separate 2D (polygon) and 3D (polyhedron) implementations
- Essential for physics collision response

**SAT (Separating Axis Theorem)**
- Tests overlap by projecting onto potential separating axes
- More accurate than GJK for polygonal contacts
- Used for: Cuboid-Cuboid, Cuboid-Triangle, Cuboid-SupportMap
- Tests face normals and edge cross products
- Deterministic (no iteration)

#### Contact Manifolds

**`ContactManifold`** - Persistent contact tracking
- Maintains multiple contact points over time
- Spatial coherence for warm-starting solvers
- Contact prediction for CCD
- Custom data storage (generic parameters)
- Essential for stable physics simulation

**`TrackedContact`** - Individual contact in manifold
- Positions on both shapes
- Penetration distance (negative = overlap)
- Feature IDs for tracking
- User data for custom physics information

**`ContactManifoldsWorkspace`** - Reusable workspace
- Avoids allocations across frames
- Caches intermediate results
- Significantly improves performance

#### Query Dispatchers

**`QueryDispatcher` trait** - Extensibility mechanism
- Allows custom shape types
- Fallback chaining support
- Override default algorithms
- Implement once, works with all queries

**`DefaultQueryDispatcher`** - Built-in implementation
- Handles all standard shape pairs
- Algorithm selection (GJK, EPA, SAT, specialized)
- Thread-safe (Send + Sync)
- Optimized for each shape combination

**`PersistentQueryDispatcher`** - With contact manifolds
- Maintains contact history
- Better physics stability
- Spatial coherence exploitation

#### Split Operations

**`SplitResult<T>`** - Shape splitting by plane
- Divides shapes into pieces
- Returns: `Pair`, `Negative`, or `Positive`
- Supported: Aabb, Segment, TriMesh
- Use cases: spatial partitioning, mesh slicing

**`IntersectResult<T>`** - Plane intersection geometry
- Computes cross-sectional geometry
- Can produce multiple components
- Use cases: contour generation, cross-sections

### 3. `bounding_volume` - Bounding Volumes

Bounding volumes for broad-phase collision detection and spatial queries.

#### Aabb - Axis-Aligned Bounding Box

**Most common bounding volume**:
- Defined by min and max corners
- Edges parallel to coordinate axes
- Extremely fast intersection tests (6 comparisons in 3D)
- Must be recomputed when objects rotate

**Key Methods**:
- Constructors: `new()`, `from_half_extents()`, `from_points()`
- Queries: `intersects()`, `contains()`, `contains_local_point()`
- Properties: `center()`, `half_extents()`, `extents()`, `volume()`
- Modifications: `merged()`, `loosened()`, `tightened()`, `transform_by()`
- Utilities: `split_at_center()`, `vertices()`, `intersection()`

**Use Cases**:
- Broad-phase collision detection
- BVH and spatial partitioning
- Frustum culling
- Quick rejection tests

#### BoundingSphere - Spherical Bounding Volume

**Rotation-invariant alternative to AABB**:
- Defined by center point and radius
- No recomputation needed when rotating
- Slightly less tight than AABB for boxes
- Simple intersection test (distance check)

**Key Methods**:
- Constructor: `new(center, radius)`
- Queries: `intersects()`, `contains()`
- Properties: `center()`, `radius()`
- Modifications: `merged()`, `loosened()`, `transform_by()`

**When to Use**:
- Objects that rotate frequently
- Spherical or near-spherical objects
- Simple broad-phase needs

#### BoundingVolume Trait

Common interface for all bounding volumes:
- `center()`, `intersects()`, `contains()`
- `merge()`, `merged()` - Combine bounding volumes
- `loosen()`, `loosened()` - Expand by margin
- `tighten()`, `tightened()` - Contract by margin

### 4. `partitioning` - Spatial Partitioning

Acceleration structures for efficient spatial queries.

#### Bvh - Bounding Volume Hierarchy

**Most important acceleration structure**:
- Binary tree with AABBs at nodes
- 10x-100x speedup for complex scenes
- Supports both static and dynamic scenes
- Replaced `Qbvh` in v0.22.0

**Construction**:
- `new_from_leaves(strategy, aabbs)` - Build from AABB collection
- **Binned strategy**: Fast build, good quality, O(n log n)
- **PLOC strategy**: Slower build, better quality
- Incremental: `insert()` one at a time

**Dynamic Updates**:
- `insert(leaf_id, aabb)` - Add or update a leaf (O(log n))
- `remove(leaf_id)` - Remove a leaf (O(log n))
- `refit()` - Update all AABBs after bulk movement (O(n), very fast)
- `optimize_incremental(budget)` - Gradually improve tree quality
- Efficient for moving objects, changing scenes

**Queries**:
- `intersect_aabb(aabb)` - Find all leaves intersecting AABB
- `cast_ray(ray)` - Ray casting through the tree
- `traverse(visitor)` - Custom traversal with early exit
- `leaves()` - Iterate all leaves
- `leaf_pairs()` - Potential collision pairs

**Performance**:
- Construction: O(n log n)
- Queries: O(log n) average with effective pruning
- Refit: O(n) but very fast (just AABB updates)
- Memory: O(n) with ~2n nodes for n leaves

**Use Cases**:
- Broad-phase collision detection
- Ray casting in complex scenes
- Frustum culling for rendering
- Spatial queries (nearby objects)

**Best Practices**:
- Use `refit()` for moving objects (not rebuild)
- Call `optimize_incremental()` every 5-10 frames
- Reuse `BvhWorkspace` across queries
- Profile to choose right construction strategy

### 5. `mass_properties` - Mass Properties

Computation of physical properties for rigid body dynamics.

#### MassProperties Struct

**Defines how objects respond to forces**:
- `local_com` - Center of mass (local coordinates)
- `inv_mass` - Inverse mass (0.0 = infinite/immovable)
- `inv_principal_inertia` - Inverse angular inertia
- `principal_inertia_local_frame` - Inertia axes (3D only)

**Why Inverse Values?**
- Infinite mass/inertia represented as zero
- Avoids division in physics loop (multiply by inverse)
- More numerically stable for heavy objects

**Constructors for All Shapes**:
```rust
// Simple shapes
MassProperties::from_ball(density, radius)
MassProperties::from_cuboid(density, half_extents)
MassProperties::from_capsule(density, a, b, radius)

// 3D shapes
MassProperties::from_cylinder(density, half_height, radius)
MassProperties::from_cone(density, half_height, radius)

// Polygonal shapes
MassProperties::from_convex_polygon(density, points) // 2D
MassProperties::from_convex_polyhedron(density, vertices, indices) // 3D

// Complex shapes
MassProperties::from_trimesh(density, vertices, indices)
MassProperties::from_compound(density, shapes)
MassProperties::from_voxels(density, voxels)
```

**Operations**:
- Add/subtract mass properties (parallel axis theorem)
- `transform_by()` - Apply rigid transformation
- `set_mass()` - Change mass with automatic inertia scaling
- `mass()`, `principal_inertia()` - Accessor methods

**Integration with Physics**:
- Directly usable in rigid body dynamics
- Center of mass for force application
- Inertia tensor for angular dynamics
- Inverse values for constraint solving

### 6. `transformation` - Mesh Transformations

Algorithms for mesh processing, analysis, and generation.

#### Convex Hull

**`convex_hull(points)`** - Compute convex hull
- 2D: Returns vertices in CCW order
- 3D: Returns vertices + triangle indices (CCW winding)
- Algorithm: Quickhull (O(n log n) average)
- Output: Minimal vertex set + connectivity

**`try_convex_hull(points)`** - Safe version
- Returns `Result` instead of panicking
- Handles degenerate inputs gracefully
- Errors: `IncompleteInput`, `MissingSupportPoint`, etc.

#### Polygon Operations

**Boolean Operations** (2D):
- `polygons_intersection()` - Intersect concave polygons
- `convex_polygons_intersection()` - Faster for convex
- Returns intersection points with connectivity
- Handles multiple connected components
- Tolerances for numerical robustness

**Convex Decomposition** (2D):
- `hertel_mehlhorn()` - Decompose into convex polygons
- Used by `Compound::decompose_trimesh()`
- Minimizes number of pieces

#### Mesh Operations

**Mesh Intersection** (3D, requires `spade` feature):
- `intersect_meshes()` - Boolean intersection of two meshes
- Requires: topology + pseudo-normals (`ORIENTED` flag)
- Returns new mesh representing intersection
- Robust to numerical issues
- Tolerances for robustness

**Mesh Generation**:
- `to_trimesh` module - Convert shapes to triangle meshes
- Quality control via subdivision parameters
- Supported for: Ball, Capsule, Cone, Cylinder, Cuboid, HeightField, etc.

#### VHACD - Convex Decomposition

**V-HACD Algorithm** - Decompose concave meshes into convex parts:
- Enables efficient collision for complex shapes
- Configurable quality vs. performance tradeoff
- Three-stage: voxelization → decomposition → hull generation

**`VHACDParameters`**:
- `concavity` - Quality vs. part count (most important)
- `resolution` - Voxel grid resolution (32-64 for games, 100-400 for quality)
- `max_convex_hulls` - Limit number of parts
- `alpha`, `beta` - Symmetry and revolution biases
- Presets: default, high-quality, fast, game-ready

**Workflow**:
1. `VHACD::decompose()` - Run decomposition
2. `compute_convex_hulls()` - Generate approximate hulls (fast)
3. OR `compute_exact_convex_hulls()` - Generate precise hulls (slower)
4. Create `Compound` from results

**Use Cases**:
- Complex rigid bodies in physics
- Destructible objects
- Character collision meshes
- Optimized static geometry

#### Voxelization

**`VoxelSet`** - Convert meshes to voxels:
- Sparse storage format
- `with_voxel_size()` - Specify voxel dimensions
- `voxelize()` - Specify grid resolution
- `FillMode`: `SurfaceOnly` vs `FloodFill` (interior)

**Use Cases**:
- Volume computation
- Collision proxies for complex shapes
- Spatial queries and sampling
- GPU-friendly representation

**`VoxelizedVolume`** - Dense storage (rarely used)

#### Transformation Utilities

**Mesh Generation Primitives** (`utils` module):
- `push_circle()` - Generate circle vertices
- `push_ring_indices()` - Connect two circles
- `push_rectangle_indices()` - Quad triangulation
- `transform()`, `transformed()` - Apply rigid transformations
- `scaled()` - Apply non-uniform scaling

**Convex Hull Utilities**:
- `support_point_id()` - Find furthest point
- `normalize()` - Center and scale point cloud

### 7. `utils` - Utility Functions

Miscellaneous geometric and mathematical utilities.

#### Geometric Utilities

**`center(points)`** - Compute centroid
- Average of all point coordinates
- O(n) linear scan

**`median(values)`** - Median of slice
- Sorts in-place (side effect!)
- Handles odd/even count correctly

**`point_in_poly2d(point, vertices)`** - 2D point containment
- Winding number algorithm
- Works for concave and self-intersecting polygons
- O(n) where n = vertex count

**`point_in_convex_poly2d(point, vertices)`** - Faster convex version
- Half-space tests
- Requires counter-clockwise winding

#### Sorting and Pairing

**`sort2()`, `sort3()`** - Small array sorting
- Minimal comparisons (3-5 for sort3)
- Returns references to sorted elements
- Useful for triangle vertex ordering

**`SortedPair<T>`** - Canonical pair ordering
- Always stores smaller element first
- Perfect for HashMap keys (edges, vertex pairs)
- Implements Hash, Eq, Ord

#### Math Utilities

**`SdpMatrix2`, `SdpMatrix3`** - Symmetric definite-positive matrices
- Used for inertia tensors
- Efficient storage and operations

**`IsometryOps`** - Isometry extensions
- Absolute transforms for AABBs
- Inverse operations

**Data Structures**:
- `VecMap` - Map optimized for small integer keys
- `hashmap`, `hashset` - Re-exports with deterministic hashing

**Morton Codes** - Z-order space-filling curves
- Spatial indexing and sorting
- Improves cache coherence

### 8. Error Handling

Parry provides clear error types with actionable solutions:

**`Unsupported`** - Query not implemented for shape pair
- Some combinations not supported
- Solutions: decompose shapes, use BVH traversal, custom dispatcher

**`ConvexHullError`** - Convex hull computation failures
- `IncompleteInput` - Too few points (<3)
- `MissingSupportPoint` - Degenerate/collinear points
- Solutions: add more points, check for duplicates

**`MeshIntersectionError`** - Mesh boolean operation failures
- `MissingTopology` - Need `TriMeshFlags::ORIENTED`
- `MissingPseudoNormals` - Enable topology first
- Solutions: set proper flags when creating mesh

**`TopologyError`** - Mesh topology validation
- `BadTriangle` - Degenerate triangle (duplicate vertices)
- `BadAdjacentTrianglesOrientation` - Inconsistent winding
- Solutions: use `DELETE_DEGENERATE_TRIANGLES` flag

**`TriMeshBuilderError`** - Mesh construction errors
- `EmptyIndices` - No triangles provided
- `TopologyError` - Wrapped topology error

**`PolygonsIntersectionError`** - Polygon operation failures
- `InfiniteLoop` - Degenerate input detected

All errors include detailed documentation explaining causes and solutions.

## Feature Flags

### Required Features (Auto-enabled)

- `dim2` / `dim3` - Compile for 2D or 3D (mutually exclusive)
- `f32` / `f64` - Floating-point precision (mutually exclusive)

### Optional Features

**Standard Library**:
- `std` (default) - Use Rust standard library
  - Enables: `alloc`, `slab`, `simba/std`, `ena`
- `alloc` - Use `alloc` crate for heap allocations
  - Required for: composite shapes, BVH, transformations, meshes

**Serialization**:
- `serde-serialize` - Serde support for all shapes
- `rkyv-serialize` - Zero-copy deserialization with `rkyv`
- `bytemuck-serialize` - `bytemuck::Pod` for simple shapes

**Performance**:
- `simd-stable` - SIMD acceleration using `wide` (stable Rust, ~4x faster)
- `simd-nightly` - SIMD using `portable_simd` (nightly Rust)
- `parallel` - Parallel BVH traversal (uses `rayon`)
- `enhanced-determinism` - Deterministic hashing and libm math
  - **Incompatible with SIMD features**

**Geometry Features**:
- `spade` (default) - Delaunay triangulation for mesh intersection
- `wavefront` (3D only) - Export meshes as OBJ files via `TriMesh::to_obj_file()`

**Advanced**:
- `improved_fixed_point_support` - Better fixed-point math support

## Type System & Generics

### The `Real` Type

Parry uses a type alias `Real` that resolves to either `f32` or `f64` based on features:

```rust
#[cfg(feature = "f32")]
pub use f32 as Real;

#[cfg(feature = "f64")]
pub use f64 as Real;
```

### Math Type Aliases

The `math` module provides dimension-agnostic type aliases:

**Common Types** (2D and 3D):
- `Point<N>` - Point2/Point3
- `Vector<N>` - Vector2/Vector3
- `Isometry<N>` - Isometry2/Isometry3 (rigid transformation)
- `Rotation<N>` - UnitComplex/UnitQuaternion
- `Matrix<N>` - Matrix2/Matrix3

**Constants**:
- `DIM` - Dimension (2 or 3)
- `DEFAULT_EPSILON` - Tolerance for geometric operations (f32::EPSILON or f64::EPSILON)

### SIMD Support

When SIMD is enabled (`simd-stable` or `simd-nightly`):
- `SimdReal` - 4-lane SIMD float (`f32x4` or `f64x4`)
- `SimdBool` - 4-lane SIMD boolean mask
- `SIMD_WIDTH` - Number of lanes (4)

Without SIMD:
- `SimdReal` - Scalar float (f32 or f64)
- `SimdBool` - Scalar bool
- `SIMD_WIDTH` - 1

## Common Usage Patterns

### Basic Distance Query

```rust
use parry3d::shape::Ball;
use parry3d::query;
use na::Isometry3;

let ball1 = Ball::new(0.5);
let ball2 = Ball::new(1.0);

let pos1 = Isometry3::identity();
let pos2 = Isometry3::translation(5.0, 0.0, 0.0);

let distance = query::distance(&pos1, &ball1, &pos2, &ball2).unwrap();
// distance = 5.0 - 0.5 - 1.0 = 3.5
```

### Ray Casting

```rust
use parry3d::shape::Cuboid;
use parry3d::query::{Ray, RayCast};
use na::{Point3, Vector3, Isometry3};

let cube = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
let ray = Ray::new(Point3::new(0.0, 0.0, -5.0), Vector3::z());

if let Some(toi) = cube.cast_ray(&Isometry3::identity(), &ray, 100.0, true) {
    println!("Hit at t = {}", toi);
    let hit_point = ray.point_at(toi);
    println!("Hit point: {:?}", hit_point);
}
```

### Contact Computation

```rust
use parry3d::query;
use parry3d::shape::{Ball, Cuboid};
use na::{Isometry3, Vector3};

let ball = Ball::new(0.5);
let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));

let pos1 = Isometry3::identity();
let pos2 = Isometry3::translation(2.0, 0.0, 0.0);

if let Ok(Some(contact)) = query::contact(&pos1, &ball, &pos2, &cuboid, 0.0) {
    println!("Contact normal: {}", contact.normal1);
    println!("Penetration: {}", -contact.dist);
}
```

### Shape Casting (Continuous Collision)

```rust
use parry3d::query::{cast_shapes, ShapeCastOptions};
use parry3d::shape::Ball;
use na::{Isometry3, Vector3};

let ball1 = Ball::new(0.5);
let ball2 = Ball::new(0.5);

let pos1 = Isometry3::identity();
let pos2 = Isometry3::translation(10.0, 0.0, 0.0);

let vel1 = Vector3::new(2.0, 0.0, 0.0); // Moving at speed 2
let vel2 = Vector3::zeros();

let options = ShapeCastOptions::default();

if let Ok(Some(hit)) = cast_shapes(&pos1, &vel1, &ball1, &pos2, &vel2, &ball2, options) {
    println!("Will collide at t = {}", hit.time_of_impact);
    // Move objects to time of impact to prevent tunneling
}
```

### BVH for Broad-Phase

```rust
use parry3d::partitioning::{Bvh, BvhBuildStrategy};
use parry3d::bounding_volume::Aabb;
use na::Point3;

// Create AABBs for 1000 objects
let mut aabbs = Vec::new();
for i in 0..1000 {
    let pos = Point3::new(i as f32, 0.0, 0.0);
    aabbs.push(Aabb::from_half_extents(pos, Vector3::new(0.5, 0.5, 0.5)));
}

// Build BVH (10x-100x faster than brute force)
let bvh = Bvh::new_from_leaves(BvhBuildStrategy::default(), aabbs);

// Find all objects in a region
let query = Aabb::new(Point3::new(40.0, -1.0, -1.0), Point3::new(60.0, 1.0, 1.0));
for leaf_id in bvh.intersect_aabb(&query) {
    println!("Object {} in region", leaf_id);
}
```

### Building a Compound Shape

```rust
use parry3d::shape::{Compound, SharedShape, Ball, Cuboid};
use na::{Isometry3, Vector3};

// Create a dumbbell (two balls connected by a bar)
let shapes = vec![
    (Isometry3::translation(-2.0, 0.0, 0.0), SharedShape::new(Ball::new(1.0))),
    (Isometry3::identity(), SharedShape::new(Cuboid::new(Vector3::new(2.0, 0.2, 0.2)))),
    (Isometry3::translation(2.0, 0.0, 0.0), SharedShape::new(Ball::new(1.0))),
];

let compound = Compound::new(shapes);
// Internal BVH automatically built for efficient queries
```

### Computing Mass Properties

```rust
use parry3d::shape::Capsule;
use na::Point3;

// Character capsule: 1.8m tall, 0.3m radius
let capsule = Capsule::new(
    Point3::new(0.0, -0.9, 0.0),
    Point3::new(0.0, 0.9, 0.0),
    0.3
);

// Human density ~985 kg/m³
let mass_props = capsule.mass_properties(985.0);

println!("Mass: {} kg", mass_props.mass());
println!("Center of mass: {:?}", mass_props.local_com);
println!("Inertia: {:?}", mass_props.principal_inertia());
```

## Dependencies

**Core Dependencies**:
- `nalgebra` (0.34) - Linear algebra (vectors, matrices, transformations)
- `simba` (0.9) - SIMD abstractions
- `num-traits` (0.2) - Numeric traits
- `approx` (0.5) - Approximate equality testing
- `either` (1) - Either type for shape variants

**Optional Dependencies**:
- `serde` (1.0) - Serialization (`serde-serialize` feature)
- `rkyv` (0.7.41) - Zero-copy serialization (`rkyv-serialize` feature)
- `bytemuck` (1) - Pod/Zeroable traits (`bytemuck-serialize` feature)
- `rayon` (1) - Data parallelism (`parallel` feature)
- `spade` (2.9) - Delaunay triangulation (`spade` feature, default)
- `obj` (0.10.2) - OBJ file export (`wavefront` feature, 3D only)
- `hashbrown` (0.16) - HashMap for `alloc` mode
- `indexmap` (2) - Ordered map (`enhanced-determinism` feature)
- `slab` (0.4) - Slab allocator (`std` feature)
- `ena` (0.14.3) - Union-find (`std` feature)

**Internal**:
- `downcast-rs` (2) - Trait object downcasting
- `bitflags` (2.3) - Bitflag enums
- `arrayvec` (0.7) - Stack-allocated vectors
- `smallvec` (1) - Small vector optimization
- `ordered-float` (5) - Ordered floating point
- `thiserror` (2) - Error derive macros
- `log` (0.4) - Logging facade
- `foldhash` (0.2) - Fast hashing
- `rstar` (0.12.0) - R*-tree (internal use)
- `glam` (0.30.4) - SIMD math (with SIMD features only)

## Performance Tips

### Shape Selection
- **Ball**: Fastest shape for all queries
- **Cuboid**: Very fast, especially with SAT
- **Capsule**: Fast, smooth collisions, great for characters
- **ConvexPolygon/Polyhedron**: Moderate (uses GJK/EPA or SAT)
- **TriMesh/Compound**: Slower (requires BVH traversal)
- **Voxels**: Moderate to slow depending on resolution

### Algorithm Selection
- **Distance only**: Use `distance()` (cheapest)
- **Overlap check**: Use `intersection_test()` (faster than contact)
- **Need normals**: Use `contact()` or `closest_points()`
- **Moving objects**: Use `cast_shapes()` for CCD

### BVH Optimization
- Build once, query many times
- Use `refit()` for moving objects (not rebuild)
- Call `optimize_incremental()` periodically
- Choose construction strategy based on use case
- Reuse `BvhWorkspace` across queries

### SIMD Usage
- Enable `simd-stable` for ~4x speedup
- Works best with batch operations
- Greatest benefit for BVH ray casting
- Incompatible with `enhanced-determinism`

### Memory Optimization
- Use `SharedShape` to share geometry
- Reuse workspaces (BVH, contact manifolds)
- Prefer `refit()` over rebuild for BVH
- Use sparse `Voxels` storage for large grids

## Testing

The codebase includes extensive tests:
- **Unit tests** in individual modules
- **Integration tests** in `crates/*/tests/`
- **Example programs** in `crates/*/examples/`
- **Doc-tests**: 300+ examples in documentation (all passing)
- Test scenarios for edge cases (penetration, degenerate geometry, etc.)

**Test Coverage Areas**:
- Distance and contact queries between all shape pairs
- Ray casting accuracy
- Convex hull computation robustness
- EPA convergence in degenerate cases
- BVH correctness after updates
- Mesh intersection robustness
- Time-of-impact accuracy
- Voxelization correctness
- Mass properties calculations

## Development Guidelines

### Code Style

- Uses `#![no_std]` compatible code where possible
- Extensive use of `#[cfg(feature = "...")]` for conditional compilation
- Clippy lints enforced: `alloc_instead_of_core`, `std_instead_of_alloc`, `std_instead_of_core`
- Documentation required (`#![warn(missing_docs)]`)
- All public items now have comprehensive documentation

### Doc-Test Feature Gates

**CRITICAL**: Because Parry compiles into four separate crates (`parry2d`, `parry3d`, `parry2d-f64`, `parry3d-f64`), doc-tests that use crate-specific imports **must** check for **both** dimension and precision features.

#### Required Pattern

Doc-tests using `parry2d::` or `parry3d::` imports MUST use:

```rust
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Ball;
/// // ... example code ...
/// # }
/// ```
```

#### Common Mistakes to Avoid

❌ **WRONG** - Only checking dimension:
```rust
/// ```rust
/// # #[cfg(feature = "dim3")]  // INCOMPLETE!
/// use parry3d::shape::Ball;
/// ```
```

❌ **WRONG** - No feature check at all:
```rust
/// ```rust
/// use parry3d::shape::Ball;  // Will fail for parry3d-f64!
/// ```
```

❌ **WRONG** - Cfg attribute outside code block:
```rust
/// ```rust
/// // ... code ...
/// # }
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {  // STRAY LINE!
///
/// ## Next Section
```

✅ **CORRECT** - Both dimension and precision:
```rust
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Ball;
/// // ... example code ...
/// # }
/// ```
```

#### Why This Matters

- **Dimension feature alone is insufficient**: Both `parry3d` (f32) and `parry3d-f64` (f64) enable the `dim3` feature
- **Doc-tests need both checks**: The test needs to know which specific crate to import from
- **Compilation will fail otherwise**: Without both checks, imports will fail when running doc-tests for f64 variants
- **Cfg must be inside code blocks**: Cfg attributes should only appear between the opening ````rust` and closing ``` fence markers

#### Feature Combinations

- `parry2d`: Requires `#[cfg(all(feature = "dim2", feature = "f32"))]`
- `parry3d`: Requires `#[cfg(all(feature = "dim3", feature = "f32"))]`
- `parry2d-f64`: Requires `#[cfg(all(feature = "dim2", feature = "f64"))]`
- `parry3d-f64`: Requires `#[cfg(all(feature = "dim3", feature = "f64"))]`

Note: Currently, all doc-tests use `f32` variants for consistency and simplicity.

### Performance Considerations

1. **SIMD**: Enable `simd-stable` or `simd-nightly` for 4x speedup on supported operations
2. **Parallelism**: Use `parallel` feature for large BVH traversals
3. **BVH Strategy**: Choose appropriate build strategy:
   - Binned: Faster build (O(n log n)), good quality
   - PLOC: Slower build, better quality
   - Incremental: For dynamic scenes with frequent changes
4. **Contact Manifolds**: Reuse `ContactManifoldsWorkspace` to avoid allocations
5. **Shape Selection**: Use simpler shapes when possible (Ball > Capsule > Cuboid > ConvexMesh > TriMesh)
6. **Spatial Coherence**: Exploit frame-to-frame coherence with persistent data structures

### Common Pitfalls

1. **Epsilon Handling**: Use `DEFAULT_EPSILON` for geometric comparisons
2. **CCD Thickness**: Ensure shapes have appropriate `ccd_thickness` for continuous collision
3. **Polygon Orientation**: 2D polygons must be counter-clockwise
4. **Mesh Quality**: Invalid meshes can cause panics; use `try_convex_hull` when unsure
5. **SIMD + Determinism**: Cannot enable both `simd-*` and `enhanced-determinism`
6. **Half-Extents**: Cuboids use half-extents, not full dimensions (common source of bugs)
7. **Prediction Distance**: Contact prediction helps CCD but can give false positives
8. **BVH Refitting**: Call `refit()` after updates, not full rebuild
9. **TriMesh Flags**: Set `ORIENTED` flag if you need topology or mesh intersection
10. **Voxel Resolution**: Higher resolution = more accurate but much slower and more memory

## Integration with Rapier

Parry is the geometric foundation for **Rapier**, a full-featured physics engine by Dimforge:
- Rapier uses Parry for all collision detection
- Parry shapes can be directly used as Rapier colliders
- Parry's `MassProperties` feed into Rapier's rigid body dynamics
- The BVH is used for broad-phase collision detection in Rapier
- Contact manifolds provide contact persistence for stable simulation
- Shape casting enables continuous collision detection

## API Complexity Guide

**Level 1 - Beginner** (Start here):
- Basic shapes: Ball, Cuboid, Capsule
- Simple queries: `distance()`, `intersection_test()`
- Ray casting with `Ray`
- Point queries with `PointQuery`

**Level 2 - Intermediate**:
- Composite shapes: Compound, TriMesh
- Contact queries: `contact()`, `ContactManifold`
- Shape casting: `cast_shapes()`
- Bounding volumes: Aabb, BoundingSphere
- Mass properties for physics

**Level 3 - Advanced**:
- BVH construction and optimization
- Custom query dispatchers
- Mesh operations (intersection, decomposition)
- Voxelization and voxel shapes
- Nonlinear shape casting
- Algorithm internals (GJK, EPA, SAT)

**Level 4 - Expert**:
- Custom shape implementation
- Support map implementation
- BVH traversal optimization
- Numerical stability tuning
- SIMD optimization

## Resources

- **Documentation**: https://parry.rs/docs
- **API Docs**: https://docs.rs/parry2d and https://docs.rs/parry3d
- **Discord**: https://discord.gg/vt9DJSW
- **Blog**: https://www.dimforge.com/blog
- **Examples**: See `crates/parry*/examples/` directory
- **Source Documentation**: Every public function has comprehensive docs with examples

## Getting Started Example

```rust
// In Cargo.toml:
// [dependencies]
// parry3d = "0.25"
// nalgebra = "0.34"

extern crate nalgebra as na;

use na::{Isometry3, Point3, Vector3};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cuboid;

fn main() {
    let cube = Cuboid::new(Vector3::new(1.0f32, 1.0, 1.0));
    let ray = Ray::new(Point3::new(0.0f32, 0.0, -1.0), Vector3::z());

    assert!(cube.intersects_ray(&Isometry3::identity(), &ray, f32::MAX));
}
```

## Documentation Quality

The Parry codebase now features:

- **World-class documentation**: Every public function documented
- **300+ doc-test examples**: All tested and passing
- **Beginner-friendly**: Concepts explained from first principles
- **Comprehensive**: Covers theory, practice, and use cases
- **Cross-referenced**: Extensive linking between related items
- **Performance notes**: Big-O complexity and optimization tips
- **Error handling**: All error types with solutions
- **Real-world examples**: Practical scenarios, not toy examples

See `DOCUMENTATION_COMPLETE.md` for full details on the documentation enhancement effort.

## Summary

Parry is a mature, high-performance collision detection library that:
- **Supports both 2D and 3D** with a shared codebase
- **Provides comprehensive shape types** (25+ including compositions)
- **Offers advanced features** like SIMD, parallelism, voxel support, mesh operations
- **Maintains API stability** while continuously improving performance
- **Serves as the foundation** for the Rapier physics engine
- **Has world-class documentation** making it accessible to all skill levels

The library is **production-ready** and **actively maintained**, with regular updates addressing performance, robustness, and new features. Whether you're building a game, robotics application, CAD software, or scientific simulation, Parry provides the collision detection primitives you need with excellent performance and a clean, well-documented API.
