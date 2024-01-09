use crate::numerics::{float::{Float, Float3}, rand::FloatRng};
use serde_derive::{Deserialize, Serialize};
use std::ops::{Deref, DerefMut, Mul};


// ===============================================================================================
// Data describing an electronic shell, and some related utilities.
// ===============================================================================================

#[cfg_attr(feature = "python", repr(C))]
#[derive(Clone, Copy, Default, Deserialize, PartialEq, Serialize)]
pub struct ElectronicShell {
    /// Shell occupancy number (might be fractionnal for mixtures).
    pub occupancy: Float,
    /// Shell ionisation energy, in MeV.
    pub energy: Float,
    /// Average momentum of shell electron(s), in MeV/_c_.
    pub momentum: Float,
}

impl ElectronicShell {
    /// Samples the momentum, in MeV/_c_, of a shell electron.
    pub fn sample<T: FloatRng>(&self, rng: &mut T) -> Float3 {
        loop {
            const SIGMA: Float = 0.5;
            let v = Float3::new(
                SIGMA * rng.normal(),
                SIGMA * rng.normal(),
                SIGMA * rng.normal(),
            );
            let r2 = v.norm2();
            let r = r2.sqrt();
            let m = (1.0 + r) * Float::exp(-2.0 * r);
            if rng.uniform01() <= m {
                return self.momentum * v;
            }
        }
    }
}

impl<'a, 'b> Mul<&'b Float> for &'a ElectronicShell {
    type Output = ElectronicShell;

    fn mul(self, rhs: &'b Float) -> Self::Output {
        Self::Output {
            occupancy: self.occupancy * rhs,
            energy: self.energy,
            momentum: self.momentum,
        }
    }
}

impl<'a, 'b> Mul<&'b ElectronicShell> for &'a Float {
    type Output = ElectronicShell;

    fn mul(self, rhs: &'b ElectronicShell) -> Self::Output {
        Self::Output {
            occupancy: self * rhs.occupancy,
            energy: rhs.energy,
            momentum: rhs.momentum,
        }
    }
}


// ===============================================================================================
// Electronic structure of a material, wrapping a collection of shells.
// ===============================================================================================

#[derive(Clone, Default, Deserialize, PartialEq, Serialize)]
pub struct ElectronicStructure (pub(crate) Vec<ElectronicShell>);

// Public API.
impl ElectronicStructure {
    /// Returns the total charge (i.e., the sum of shells occupancies).
    pub fn charge(&self) -> Float {
        self
            .iter()
            .map(|shell: &ElectronicShell| shell.occupancy)
            .sum()
    }

    /// Creates a new structure from the weighted sum of others. Note that atomic weights are
    /// expected, not mass fractions.
    ///
    /// # Examples
    ///
    /// ```
    /// # use goupil::physics::ElectronicStructure;
    /// # let hydrogen = ElectronicStructure::default();
    /// # let oxygen = ElectronicStructure::default();
    /// #
    /// let water = ElectronicStructure::from_others(&[
    ///     (2.0, &hydrogen),
    ///     (1.0, &oxygen),
    /// ]);
    /// ```
    pub fn from_others(composition: &[(Float, &Self)]) -> Self {
        let mut electrons = ElectronicStructure::default();
        for (weight, structure) in composition.iter() {
            for shell in structure.iter() {
                electrons.push(weight * shell);
            }
        }
        electrons
    }

    /// Samples a shell according to occupancy numbers.
    pub fn sample<T: FloatRng>(&self, rng: &mut T) -> Option<&ElectronicShell> {
        let charge = self.charge();
        if charge <= 0.0 { return None }
        let u = charge * rng.uniform01();
        let mut z = 0.0;
        for shell in self.iter() {
            z += shell.occupancy;
            if u <= z { return Some(shell) }
        }
        if (u <= charge) && (self.len() > 0) {
            return self.last();
        } else {
            return None;
        }
    }
}

impl Deref for ElectronicStructure {
    type Target = Vec<ElectronicShell>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for ElectronicStructure {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
