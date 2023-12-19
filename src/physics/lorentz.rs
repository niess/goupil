use crate::numerics::{float::{Float, Float3}, sq};
use std::ops::Neg;


// ===============================================================================================
// Four momentum data, with a few utilities.
//
// Four momentum is defined using invariant mass instead of energy. This is to avoid numeric issues
// in the evaluation of M^2 = E^2 - P^2, where P is large. It might not be the most efficient CPU
// wise, but it is numerically robust.
// ===============================================================================================
#[derive(Copy, Clone, PartialEq)]
pub struct FourMomentum {
    pub mass: Float,
    pub momentum: Float3,
}

impl FourMomentum {
    pub const fn new(mass:Float, momentum: Float3) -> Self {
        Self {mass, momentum}
    }

    pub fn energy2(self) -> Float {
        sq!(self.mass) + self.momentum.norm2()
    }

    pub fn energy(self) -> Float {
        if self.mass == 0.0 { self.momentum.norm() } else { self.energy2().sqrt() } 
    }
}


// ===============================================================================================
// Lorentz transformation.
//
// Tranform to / from the rest frame of a given particle (expressed by its four momentum). Note
// that CM frames of composite systems are not handled.
// ===============================================================================================
#[derive(Copy, Clone, PartialEq)]
pub struct LorentzTransformation {
    beta: Float3,
    gamma: Float,
}

impl LorentzTransformation {
    pub fn new(four_momentum: FourMomentum) -> Self {
        let energy = four_momentum.energy();
        let gamma = energy / four_momentum.mass;
        let beta = four_momentum.momentum / energy;
        Self {beta, gamma}
    }

    // In-place inversion.
    pub fn reverse(&mut self) { self.beta.reverse(); }

    // Apply Lorentz tranformation.
    pub fn transform(self, four_momentum: FourMomentum) -> FourMomentum {
        let dp = self.beta.dot(four_momentum.momentum);
        let energy = four_momentum.energy();
        let tmp = self.gamma * (self.gamma / (self.gamma + 1.0) * dp - energy);
        FourMomentum {
            mass: four_momentum.mass,
            momentum: four_momentum.momentum + tmp * self.beta
        }
    }
}

impl Neg for LorentzTransformation {
    type Output = Self;

    fn neg(self) -> Self {
        Self {beta: -self.beta, gamma: self.gamma}
    }
}


// ===============================================================================================
// Unit tests.
// ===============================================================================================
#[cfg(test)]
mod tests {
    use crate::numerics::{consts::SQRT_2, tests::assert_float_eq};
    use crate::physics::consts::ELECTRON_MASS;
    use super::*;

    #[test]
    fn four_momentum() {
        let p = FourMomentum::new(ELECTRON_MASS, Float3::new(0.0, 0.0, 2.0 * ELECTRON_MASS));
        assert_eq!(p.mass, ELECTRON_MASS);
        let tmp = 5.0 * sq!(ELECTRON_MASS);
        assert_eq!(p.energy2(), tmp);
        assert_eq!(p.energy(), tmp.sqrt());
    }

    #[test]
    fn lorentz_transformation() {
        let p0 = FourMomentum::new(ELECTRON_MASS, Float3::new(0.0, 0.0, ELECTRON_MASS));
        let t0 = LorentzTransformation::new(p0);
        assert_float_eq!(t0.gamma, SQRT_2);

        let mut t1 = -t0;
        assert!(t1.beta == -t0.beta);

        t1.reverse();
        assert!(t1 == t0);

        let p1 = t0.transform(p0);
        assert_float_eq!(p1.momentum.norm2(), 0.0);

        t1.reverse();
        let p2 = t1.transform(p1);
        assert!(p0 == p2);
    }
}
