use super::float::Float;

// Floating point constants.
#[cfg(not(feature = "f32"))]
use std::f64::consts;
#[cfg(feature = "f32")]
use std::f32::consts;

macro_rules! export_float_consts {
    ( $( $symbol:ident ),* ) => {
        $(
            pub const $symbol: Float = consts::$symbol;
        )*
    }
}

export_float_consts!(
    E,
    FRAC_1_PI,
    FRAC_1_SQRT_2,
    FRAC_2_PI,
    FRAC_2_SQRT_PI,
    FRAC_PI_2,
    FRAC_PI_3,
    FRAC_PI_4,
    FRAC_PI_6,
    FRAC_PI_8,
    LN_2,
    LN_10,
    LOG2_10,
    LOG2_E,
    LOG10_2,
    LOG10_E,
    PI,
    SQRT_2,
    TAU
);
