## Repository architecture

The architecture of this repository is a bit unusual because we are using some tricks to have both
the 2D and 3D version of Parry share the same code-base. Here are the main folders:
- **`build/`**: contains one folder per Parry crate (for the 2D, 3D, `f32`, and `f64` versions). Each
  crate has its own `Cargo.toml` file that adjusts some cargo features, and reference the `src` folder.
- **`src/`**: contains the actual `.rs` source code of the Parry geometric library.
