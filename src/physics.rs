pub(crate) mod elements;
pub(crate) mod materials;
pub(crate) mod process;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub mod consts;

pub use self::elements::AtomicElement;
pub use self::materials::{
    ElectronicShell,
    ElectronicStructure,
    MaterialDefinition,
    MaterialRecord,
    MaterialRegistry,
    WeightedElement,
    WeightedMaterial,
};
pub use self::process::{
    AbsorptionMode,
    AbsorptionCrossSection,
    ComptonCrossSection,
    ComptonCDF,
    ComptonInverseCDF,
    ComptonMethod,
    ComptonMode,
    ComptonModel,
    ComptonSampler,
    RayleighCrossSection,
    RayleighFormFactor,
    RayleighMode,
};
