pub(crate) mod float;
pub(crate) mod grids;
pub(crate) mod interpolate;
pub(crate) mod integrate;
pub(crate) mod rand;
pub(crate) mod table;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub mod consts;

pub use self::float::{Float, Float3, Float3x3};
pub use self::rand::FloatRng;


// ===============================================================================================
// Unit tests.
// ===============================================================================================
#[cfg(test)]
pub mod tests {
    macro_rules! assert_float_eq {
        ($lhs:expr, $rhs:expr $(,$epsilon:expr)?) => {
            #[allow(unused_assignments)]
            #[allow(unused_mut)]
            let mut epsilon = Float::EPSILON;
            $( epsilon = $epsilon; )?
            if ($lhs - $rhs).abs() > epsilon { assert_eq!($lhs, $rhs) }
        }
    }
    pub(crate) use assert_float_eq;

    macro_rules! assert_float3_eq {
        ($lhs:expr, $rhs:expr) => {
            assert!(($lhs.0 - $rhs.0).abs() <= Float::EPSILON);
            assert!(($lhs.1 - $rhs.1).abs() <= Float::EPSILON);
            assert!(($lhs.2 - $rhs.2).abs() <= Float::EPSILON);
        }
    }
    pub(crate) use assert_float3_eq;
}
