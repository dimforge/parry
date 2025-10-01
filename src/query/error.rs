use core::fmt;

/// Error indicating that a query is not supported between certain shapes
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Unsupported;

impl fmt::Display for Unsupported {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad("query not supported between these shapes")
    }
}

#[cfg(feature = "alloc")]
impl core::error::Error for Unsupported {}
