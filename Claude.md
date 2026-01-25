# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

Parry is a 2D and 3D geometric and collision detection library written in Rust by Dimforge. It uses an unusual multi-crate architecture to share a single codebase between 2D, 3D, f32, and f64 variants.

## Critical Architecture Pattern

**This repository uses a unique shared-source architecture:**

- **`src/`** - Contains all `.rs` source code (shared by all crates)
- **`crates/`** - Contains four separate crates, each with its own `Cargo.toml`:
  - `parry2d` - 2D with f32 precision
  - `parry3d` - 3D with f32 precision
  - `parry2d-f64` - 2D with f64 precision
  - `parry3d-f64` - 3D with f64 precision

Each crate's `Cargo.toml` points to `../../src/lib.rs` and enables specific features (`dim2`/`dim3`, `f32`/`f64`) to compile the appropriate code from the shared source.

## Common Development Commands

### Building

```bash
# Build all workspace crates
cargo build

# Build a specific crate
cargo build -p parry3d
cargo build -p parry2d

# Build with SIMD optimizations
cd crates/parry3d && cargo build --features simd-stable

# Build with all serialization features
cargo build --features bytemuck-serialize,serde-serialize,rkyv

# Build with enhanced determinism (incompatible with SIMD)
cargo build --features enhanced-determinism

# Build for WASM
cd crates/parry3d && cargo build --target wasm32-unknown-unknown
```

### Testing

```bash
# Run all tests
cargo test

# Run tests with specific features
cargo test --features wavefront
cargo test --features parallel

# Run a single test
cargo test test_name

# Run tests for a specific crate
cargo test -p parry3d

# Run doc-tests (IMPORTANT: these test all four crate variants)
cargo test --doc
```

### Code Quality

```bash
# Check formatting (required before committing)
cargo fmt -- --check

# Format code
cargo fmt

# Run clippy (warnings treated as errors in CI)
cargo clippy

# Build documentation
cargo doc --workspace --no-deps

# Build documentation with private items
cargo doc --workspace --no-deps --document-private-items
```

### Examples

Examples are located in `crates/parry2d/examples/` and `crates/parry3d/examples/`:

```bash
# Run a 3D example
cargo run -p parry3d --example getting_started

# Run a 2D example
cargo run -p parry2d --example getting_started

# List all examples
ls crates/parry3d/examples/
```

## Critical Development Constraints

### 1. Doc-Test Feature Gates

**CRITICAL**: Because the codebase compiles into four separate crates, doc-tests that import crate-specific types (like `parry3d::` or `parry2d::`) **MUST** check for BOTH dimension AND precision features:

```rust
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Ball;
/// // ... example code ...
/// # }
/// ```
```

**Common mistakes:**
- ❌ Only checking dimension: `#[cfg(feature = "dim3")]` (WRONG - will fail for f64 variants)
- ❌ No feature check at all (WRONG - will fail for different variants)
- ✅ Correct: `#[cfg(all(feature = "dim3", feature = "f32"))]`

Feature combinations:
- `parry2d`: `#[cfg(all(feature = "dim2", feature = "f32"))]`
- `parry3d`: `#[cfg(all(feature = "dim3", feature = "f32"))]`
- `parry2d-f64`: `#[cfg(all(feature = "dim2", feature = "f64"))]`
- `parry3d-f64`: `#[cfg(all(feature = "dim3", feature = "f64"))]`

### 2. Conditional Compilation Patterns

Code in `src/` uses extensive `#[cfg(feature = "...")]` attributes:

```rust
#[cfg(feature = "dim2")]
// 2D-specific code

#[cfg(feature = "dim3")]
// 3D-specific code

#[cfg(feature = "f32")]
pub use f32 as Real;

#[cfg(feature = "f64")]
pub use f64 as Real;
```

When adding new code, consider which feature combinations it applies to.

### 3. no_std Compatibility

The codebase maintains `#![no_std]` compatibility:
- Use `core::` instead of `std::` when possible
- Use `alloc::` for heap allocations (behind `alloc` feature)
- Clippy lints enforce: `alloc_instead_of_core`, `std_instead_of_alloc`, `std_instead_of_core`

### 4. SIMD vs. Determinism

**IMPORTANT**: `simd-stable`/`simd-nightly` features are incompatible with `enhanced-determinism`. Don't enable both.

## Source Code Organization

```
src/
├── lib.rs                    # Crate root, exports all modules
├── bounding_volume/          # AABB, BoundingSphere, etc.
├── mass_properties/          # Mass, inertia calculations
├── partitioning/             # BVH spatial acceleration
├── query/                    # Collision queries (distance, contact, ray-cast)
├── shape/                    # All geometric shapes (Ball, Cuboid, TriMesh, etc.)
├── transformation/           # Mesh operations, convex hull, VHACD
└── utils/                    # Utilities, math helpers
```

## Adding New Features

When adding new public APIs:

1. **Documentation**: Add comprehensive doc comments with examples
2. **Doc-tests**: Include working examples with proper feature gates
3. **Unit tests**: Add tests in the same file or `#[cfg(test)]` module
4. **Integration tests**: If needed, add to `crates/parry*/tests/`
5. **Examples**: Consider adding example to `crates/parry*/examples/`
6. **Feature compatibility**: Consider if it works with `no_std`, requires `alloc`, etc.

## Running CI Checks Locally

Before submitting a PR, ensure these pass (mirrors CI):

```bash
# Format check
cargo fmt -- --check

# Clippy
RUSTFLAGS="-D warnings" cargo clippy

# Build all crates
cargo build -p parry2d
cargo build -p parry3d

# Tests with various features
cargo test --features wavefront
cargo test --features parallel

# Documentation
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps

# WASM build
cd crates/parry3d && cargo build --target wasm32-unknown-unknown
```

## Publishing (Maintainers Only)

The `publish.sh` script publishes all four crates to crates.io:
- Creates temp directory with flattened structure
- Publishes in order: parry2d → parry2d-f64 → parry3d → parry3d-f64

## Key Design Patterns

### Type Aliases for Dimension-Agnostic Code

The `math` module provides type aliases that resolve based on features:
- `Point<N>` → `Point2` or `Point3`
- `Vector<N>` → `Vector2` or `Vector3`
- `Isometry<N>` → `Isometry2` or `Isometry3`
- `Real` → `f32` or `f64`

This allows writing dimension-agnostic algorithms.

### Shape Type Erasure

`SharedShape` wraps shapes in `Arc<dyn Shape>` for:
- Heterogeneous collections
- Cheap cloning via reference counting
- Copy-on-write with `make_mut()`

### Query Dispatchers

The `QueryDispatcher` trait provides extensibility:
- `DefaultQueryDispatcher` handles all built-in shape pairs
- Custom dispatchers can override specific combinations
- Algorithm selection (GJK, EPA, SAT, specialized) per shape pair

## Common Pitfalls

1. **Half-extents vs full dimensions**: Cuboids use half-extents (full_width = 2 * half_extent)
2. **Counter-clockwise winding**: 2D polygons must be CCW
3. **TriMesh flags**: Set `ORIENTED` flag for topology/mesh intersection features
4. **Prediction distance**: Contact prediction can give false positives if too large
5. **BVH updates**: Use `refit()` for moving objects, not full rebuild

## Resources

- **User Guide**: https://parry.rs/docs
- **API Docs**: https://docs.rs/parry3d (3D) and https://docs.rs/parry2d (2D)
- **Discord**: https://discord.gg/vt9DJSW
- **Blog**: https://www.dimforge.com/blog

For detailed architectural information, algorithm descriptions, and comprehensive API documentation, see `Claude.md` in the repository root.
