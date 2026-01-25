#!/bin/bash
# Run all parry 2D and 3D examples that use kiss3d for visualization.
# Close each window to proceed to the next example.

set -e

cd "$(dirname "$0")"

# 3D examples using kiss3d
EXAMPLES_3D=(
    aabb3d
    bounding_sphere3d
    convex_hull3d
    plane_intersection
    project_point3d
)

# 3D examples requiring the wavefront feature
EXAMPLES_3D_WAVEFRONT=(
    convex_decomposition
)

# 2D examples using kiss3d
EXAMPLES_2D=(
    aabb2d
    bounding_sphere2d
    convex_hull2d
    project_point2d
    point_in_poly2d
    polygons_intersection2d
    raycasts_animated
)

TOTAL=$((${#EXAMPLES_3D[@]} + ${#EXAMPLES_3D_WAVEFRONT[@]} + ${#EXAMPLES_2D[@]}))
echo "Running $TOTAL parry examples with kiss3d visualization..."
echo "Close each window to proceed to the next example."
echo ""

echo "=== Running 3D examples ==="
for example in "${EXAMPLES_3D[@]}"; do
    echo "=== Running: $example (parry3d) ==="
    cargo run --release -p parry3d --example "$example"
    echo ""
done

echo "=== Running 3D examples (with wavefront feature) ==="
for example in "${EXAMPLES_3D_WAVEFRONT[@]}"; do
    echo "=== Running: $example (parry3d) ==="
    cargo run --release -p parry3d --features wavefront --example "$example"
    echo ""
done

echo "=== Running 2D examples ==="
for example in "${EXAMPLES_2D[@]}"; do
    echo "=== Running: $example (parry2d) ==="
    cargo run --release -p parry2d --example "$example"
    echo ""
done

echo "All examples completed!"
