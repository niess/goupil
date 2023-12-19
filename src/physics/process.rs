pub(crate) mod absorption;
pub(crate) mod compton;
pub(crate) mod rayleigh;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::{
    absorption::{
        AbsorptionMode,
        AbsorptionCrossSection,
    },
    compton::{
        ComptonCrossSection,
        ComptonCDF,
        ComptonInverseCDF,
        ComptonMethod,
        ComptonMode,
        ComptonModel,
        ComptonSampler,
    },
    rayleigh::{
        RayleighCrossSection,
        RayleighFormFactor,
        RayleighMode,
    },
};
