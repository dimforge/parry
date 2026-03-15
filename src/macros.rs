/// Inner implementation of [`for_each_dim!`]. A tt-muncher that replaces every
/// occurrence of the sentinel `[ii]` in the token stream with `.$f` (i.e. `.x`,
/// `.y`, or `.z`) for the current dimension.
#[doc(hidden)]
macro_rules! for_each_dim_inner {
    // Base case: no more input tokens, emit accumulated output.
    ($f:ident [$($out:tt)*]) => { $($out)* };

    // Match the [ii] sentinel and replace with .field_name.
    ($f:ident [$($out:tt)*] [ii] $($rest:tt)*) => {
        for_each_dim_inner!($f [$($out)* . $f] $($rest)*)
    };

    // Recurse into `{ }` delimited groups.
    ($f:ident [$($out:tt)*] { $($inner:tt)* } $($rest:tt)*) => {
        for_each_dim_inner!($f [$($out)* { for_each_dim_inner!($f [] $($inner)*) }] $($rest)*)
    };

    // Recurse into `( )` delimited groups.
    ($f:ident [$($out:tt)*] ( $($inner:tt)* ) $($rest:tt)*) => {
        for_each_dim_inner!($f [$($out)* ( for_each_dim_inner!($f [] $($inner)*) )] $($rest)*)
    };

    // Recurse into `[ ]` delimited groups (general case, NOT [ii]).
    ($f:ident [$($out:tt)*] [ $($inner:tt)* ] $($rest:tt)*) => {
        for_each_dim_inner!($f [$($out)* [ for_each_dim_inner!($f [] $($inner)*) ]] $($rest)*)
    };

    // Pass through any other token unchanged.
    ($f:ident [$($out:tt)*] $other:tt $($rest:tt)*) => {
        for_each_dim_inner!($f [$($out)* $other] $($rest)*)
    };
}

/// Unrolls a loop over spatial dimensions (2 or 3 depending on features),
/// replacing occurrences of the sentinel `[ii]` with direct field access
/// (`.x`, `.y`, `.z`) on glam vectors. This avoids bracket indexing which
/// is not supported on SPIR-V targets.
///
/// # Usage
///
/// Basic form — every `[ii]` in the body is replaced with `.x`, `.y`, `.z`
/// in successive unrolled iterations:
/// ```ignore
/// for_each_dim! {
///     min[ii] = a[ii].min(b[ii]);
/// }
/// // Expands to:
/// //   min.x = a.x.min(b.x);
/// //   min.y = a.y.min(b.y);
/// //   min.z = a.z.min(b.z);  // only in 3D
/// ```
///
/// With a dimension index (provided as a `const usize`):
/// ```ignore
/// for_each_dim! { i =>
///     if mask & (1 << i) != 0 {
///         result[ii] = -result[ii];
///     }
/// }
/// ```
///
/// # Limitations
///
/// - The sentinel `[ii]` is replaced **everywhere** it appears in the body,
///   so avoid using `[ii]` for other purposes inside the macro.
macro_rules! for_each_dim {
    ($i:ident => $($body:tt)*) => {{
        #[cfg(feature = "dim2")]
        {
            { const $i: usize = 0; for_each_dim_inner!(x [] $($body)*); }
            { const $i: usize = 1; for_each_dim_inner!(y [] $($body)*); }
        }
        #[cfg(feature = "dim3")]
        {
            { const $i: usize = 0; for_each_dim_inner!(x [] $($body)*); }
            { const $i: usize = 1; for_each_dim_inner!(y [] $($body)*); }
            { const $i: usize = 2; for_each_dim_inner!(z [] $($body)*); }
        }
    }};
    ($($body:tt)*) => {{
        #[cfg(feature = "dim2")]
        {
            for_each_dim_inner!(x [] $($body)*);
            for_each_dim_inner!(y [] $($body)*);
        }
        #[cfg(feature = "dim3")]
        {
            for_each_dim_inner!(x [] $($body)*);
            for_each_dim_inner!(y [] $($body)*);
            for_each_dim_inner!(z [] $($body)*);
        }
    }};
}
