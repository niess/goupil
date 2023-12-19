use crate::physics::process::{
    absorption::table::AbsorptionCrossSection,
    compton::table::ComptonTable,
    rayleigh::table::RayleighTable,
};
use serde_derive::{Deserialize, Serialize};


// ===============================================================================================
// Physics table relative to a material.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct MaterialTable
{
    // Photo-absorption cross-section.
    pub(crate) absorption: Option<AbsorptionCrossSection>,

    // Compton process related tabulations.
    pub(crate) compton: ComptonTable,

    // Rayleigh process related tabulations.
    pub(crate) rayleigh: RayleighTable,
}
