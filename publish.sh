#! /bin/bash
#
# Publishes all four parry crates (parry2d, parry3d, parry2d-f64, parry3d-f64)
# to crates.io using a single `cargo publish --workspace` invocation.
#
# Why this script exists
# ----------------------
# All four crates share a single source tree at `<repo>/src`, referenced from
# each crate's manifest as `path = "../../src/lib.rs"`. That path points outside
# the crate directory, which `cargo publish` refuses to package.
#
# To work around it *only during publishing*, this script temporarily:
#   1. rewrites each crate's `[lib] path` from `../../src/lib.rs` to `src/lib.rs`, and
#   2. creates a `src` symlink in each crate pointing at the shared `../../src`.
#
# Cargo follows the symlink and bundles the real source into each `.crate`, so
# native workspace publishing works. A trap restores the manifests and removes
# the symlinks on exit (including on error or Ctrl-C), leaving the tree as it was.
#
# Extra arguments are forwarded to `cargo publish`, e.g.:
#   ./publish.sh --dry-run
#   ./publish.sh --token "$CARGO_TOKEN"
#
# Requires cargo >= 1.90 (for `cargo publish --workspace`).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

CRATES=(parry2d parry3d parry2d-f64 parry3d-f64)

# Refuse to run on a dirty tree: the only diff during publishing must be our own
# temporary edits, so the restore at the end is guaranteed to be correct.
if [ -n "$(git status --porcelain)" ]; then
    echo "error: working tree is not clean. Commit or stash changes before publishing." >&2
    exit 1
fi

backup_dir="$(mktemp -d)"

cleanup() {
    for crate in "${CRATES[@]}"; do
        # Remove the symlink we created (only if it is in fact a symlink).
        if [ -L "crates/$crate/src" ]; then
            rm -f "crates/$crate/src"
        fi
        # Restore the original manifest.
        if [ -f "$backup_dir/$crate.Cargo.toml" ]; then
            cp "$backup_dir/$crate.Cargo.toml" "crates/$crate/Cargo.toml"
        fi
    done
    rm -rf "$backup_dir"
}
trap cleanup EXIT INT TERM

# Apply the temporary symlink layout.
for crate in "${CRATES[@]}"; do
    manifest="crates/$crate/Cargo.toml"
    cp "$manifest" "$backup_dir/$crate.Cargo.toml"

    tmp="$(mktemp)"
    sed 's#path = "\.\./\.\./src/lib.rs"#path = "src/lib.rs"#' "$manifest" > "$tmp"
    mv "$tmp" "$manifest"

    ln -s ../../src "crates/$crate/src"
done

# Publish the whole workspace. `--allow-dirty` is required because our temporary
# edits make the tree dirty; the clean-tree check above keeps that safe.
cargo publish --workspace --allow-dirty "$@"
