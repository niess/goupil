use crate::numerics::{
    float::Float,
    grids::UnstructuredGrid,
    interpolate::CubicInterpolator,
    table::{Data1D, Table1D},
};
use crate::physics::consts::BARN;
use serde_derive::{Deserialize, Serialize};


// ===============================================================================================
// Photo-absorption cross-section table, and related utilities.
// ===============================================================================================

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct AbsorptionCrossSection {
    pub(crate) energies: UnstructuredGrid,
    pub(crate) values: Vec<Float>,
    pub(super) interpolator: CubicInterpolator,
}

// Public API.
impl AbsorptionCrossSection {
    pub fn energies(&self) -> &[Float] {
        self.energies.as_ref()
    }

    pub fn interpolate(&self, energy: Float) -> Float {
        self.interpolator
            .interpolate(&self.energies, &self.values, energy)
            .unwrap_or(0.0)
    }

    pub fn len(&self) -> usize {
        self.energies.len()
    }

    pub fn values(&self) -> &[Float] {
        &self.values.as_ref()
    }
}

// Private API.
impl AbsorptionCrossSection {
    // Merge multiple absorption data into a single one, taking care of preserving spectral lines.
    // Note that mole weights are expected.
    pub(crate) fn from_others(tables: &[(Float, &Self)]) -> Option<Self> {
        Table1D::merge(tables)
            .map(|(energies, values)| Self::new(energies, values))
    }

    pub(crate) fn new(energies: Vec<Float>, values: Vec<Float>) -> Self {
        let mut interpolator = CubicInterpolator::new(energies.len());
        let energies = UnstructuredGrid::from(energies);
        interpolator.initialise(&energies, &values, false);
        Self {energies, values, interpolator}
    }
}

impl Table1D for AbsorptionCrossSection {
    #[inline]
    fn interpolate(&self, x: Float) -> Float {
        Self::interpolate(self, x)
    }

    #[inline]
    fn len(&self) -> usize {
        Self::len(self)
    }

    #[inline]
    fn x(&self, i: usize) -> Float {
        self.energies[i]
    }

    #[inline]
    fn y(&self, i: usize) -> Float {
        self.values[i]
    }
}

impl From<Data1D> for AbsorptionCrossSection {
    fn from(mut data: Data1D) -> Self {
        for y in data.y.iter_mut() {
            *y *= BARN;
        }
        Self::new(data.x, data.y)
    }
}
