use crate::numerics::float::Float;
use crate::numerics::grids::UnstructuredGrid;
use crate::numerics::interpolate::CubicInterpolator;
use crate::numerics::table::{Data1D, Table1D};
use crate::physics::consts::BARN;
use serde_derive::{Deserialize, Serialize};


// ===============================================================================================
// Rayleigh process tabulations.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub(crate) struct RayleighTable {
    pub cross_section: Option<RayleighCrossSection>,
    pub form_factor: Option<RayleighFormFactor>,
}


// ===============================================================================================
// Rayleigh scattering cross-section.
// ===============================================================================================

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct RayleighCrossSection {
    pub(crate) energies: UnstructuredGrid,
    pub(crate) values: Vec<Float>,
    pub(crate) interpolator: CubicInterpolator,
}

// Public API.
impl RayleighCrossSection {
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
impl RayleighCrossSection {
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

impl Table1D for RayleighCrossSection {
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

impl From<Data1D> for RayleighCrossSection {
    fn from(mut data: Data1D) -> Self {
        for y in data.y.iter_mut() {
            *y *= BARN;
        }
        Self::new(data.x, data.y)
    }
}


// ===============================================================================================
// Rayleigh scattering form-factor, and related utilities.
// ===============================================================================================

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct RayleighFormFactor {
    pub(crate) momenta: UnstructuredGrid,
    pub(crate) values: Vec<Float>,
    pub(crate) interpolator: CubicInterpolator,
    // Enveloppe parameter(s).
    pub(crate) scale: Float,
}

// Public API.
impl RayleighFormFactor {
    pub fn interpolate(&self, momentum: Float) -> Float {
        self.interpolator
            .interpolate(&self.momenta, &self.values, momentum)
            .unwrap_or(0.0)
    }

    pub fn momenta(&self) -> &[Float] {
        self.momenta.as_ref()
    }

    pub fn len(&self) -> usize {
        self.momenta.len()
    }

    pub fn values(&self) -> &[Float] {
        &self.values.as_ref()
    }
}

// Private API.
impl RayleighFormFactor {
    pub(crate) fn enveloppe(&self, squared_momentum: Float) -> Float {
        self.values[0] * self.scale / (self.scale + squared_momentum)
    }

    pub(crate) fn from_others(tables: &[(Float, &Self)]) -> Option<Self> {
        Table1D::merge(tables)
            .map(|(energies, values)| Self::new(energies, values))
    }

    // Note: this is the inverse CDF of the squared enveloppe, as function of squared momentum.
    // That is, it returns the squared momentum for a given CDF value.
    pub(crate) fn inverse_cdf(&self, cdf: Float, q2max: Float) -> Float {
        if cdf <= 0.0 {
            0.0
        } else if cdf >= 1.0 {
            q2max
        } else {
            let a = cdf * q2max / (self.scale + q2max);
            let x = a / (1.0 - a);
            self.scale * x
        }
    }

    pub(crate) fn new(momenta: Vec<Float>, values: Vec<Float>) -> Self {
        let scale = {
            let n = momenta.len();
            let mut scale = 0.0;
            for i in 1..n {
                let x = momenta[i] * momenta[i];
                let f = values[i] / values[0];
                if f < 1.0 {
                    let b = x * f / (1.0 - f);
                    if b > scale {
                        scale = b;
                    }
                }
            }
            scale
        };
        let mut interpolator = CubicInterpolator::new(momenta.len());
        let momenta = UnstructuredGrid::from(momenta);
        interpolator.initialise(&momenta, &values, false);
        Self {momenta, values, interpolator, scale}
    }
}

impl Table1D for RayleighFormFactor {
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
        self.momenta[i]
    }

    #[inline]
    fn y(&self, i: usize) -> Float {
        self.values[i]
    }
}

impl From<Data1D> for RayleighFormFactor {
    fn from(data: Data1D) -> Self {
        Self::new(data.x, data.y)
    }
}
