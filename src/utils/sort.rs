#[cfg(feature = "dim3")]
#[inline]
/// Sorts a set of three values in increasing order.
pub fn sort2<T: PartialOrd + Copy>(a: T, b: T) -> (T, T) {
    if a > b {
        (b, a)
    } else {
        (a, b)
    }
}

/// Sorts a set of three values in increasing order.
#[inline]
pub fn sort3<'a, T: PartialOrd + Copy>(a: &'a T, b: &'a T, c: &'a T) -> (&'a T, &'a T, &'a T) {
    let a_b = *a > *b;
    let a_c = *a > *c;
    let b_c = *b > *c;

    let sa;
    let sb;
    let sc;

    // Sort the three values.
    if a_b {
        // a > b
        if a_c {
            // a > c
            sc = a;

            if b_c {
                // b > c
                sa = c;
                sb = b;
            } else {
                // b <= c
                sa = b;
                sb = c;
            }
        } else {
            // a <= c
            sa = b;
            sb = a;
            sc = c;
        }
    } else {
        // a < b
        if !a_c {
            // a <= c
            sa = a;

            if b_c {
                // b > c
                sb = c;
                sc = b;
            } else {
                sb = b;
                sc = c;
            }
        } else {
            // a > c
            sa = c;
            sb = a;
            sc = b;
        }
    }

    (sa, sb, sc)
}
