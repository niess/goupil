use anyhow::{anyhow, Result};
use crate::numerics::float::{Float, Float3};
use std::any::type_name;


// ===============================================================================================
// Density models and related routines.
// ===============================================================================================

#[derive(Copy, Clone)]
pub enum DensityModel {

    // Exponential gradient, in g / cm^3..
    Gradient {
        rho0: Float,
        origin: Float3,
        lambda: Float,
        direction: Float3,
    },

    // Uniform density, in g / cm^3.
    Uniform (Float),
}

impl DensityModel {
    // Column depth from geometrical range.
    pub fn column_depth(&self, position: Float3, direction: Float3, range: Float) -> Float {
        match self {
            Self::Gradient{rho0, origin: r0, lambda, direction: u0} => {
                let rho = rho0 * (u0.dot(position - *r0) / lambda).exp();
                let costheta = u0.dot(direction);
                if costheta.abs() < 1e-7 {
                    range * rho
                } else {
                    let lambda = lambda / costheta;
                    lambda * rho * ((range / lambda).exp() - 1.0)
                }
            },
            Self::Uniform(rho) => { *rho * range },
        }
    }

    pub fn gradient(rho0: Float, origin: Float3, lambda: Float, direction: Float3)
        -> Result<Self> {

        if rho0 <= 0.0 {
            Err(anyhow!(
                "{}: bad density (expected a strictly positive value, found {})",
                type_name::<Self>(),
                rho0,
            ))
        } else {
            Ok(Self::Gradient {rho0, origin, lambda, direction})
        }
    }

    // Geometrical range, in cm, for given column depth, in g / cm^2.
    pub fn range(&self, position: Float3, direction: Float3, column_depth: Float) -> Float {
        match self {
            Self::Gradient{rho0, origin: r0, lambda, direction: u0} => {
                let rho = rho0 * (u0.dot(position - *r0) / lambda).exp();
                if rho <= 0.0 {
                    Float::INFINITY
                } else {
                    let costheta = u0.dot(direction);
                    if costheta.abs() < 1e-7 {
                        column_depth / rho
                    } else {
                        let lambda = lambda / costheta;
                        let r = column_depth / (rho * lambda);
                        if r <= -1.0 {
                            Float::INFINITY
                        } else {
                            lambda * (1.0 + r).ln()
                        }
                    }
                }
            },
            Self::Uniform(rho) => if *rho <= 0.0 { Float::INFINITY } else { column_depth / rho },
        }
    }

    pub fn uniform(density: Float) -> Result<Self> {
        Ok(Self::Uniform (density))
    }

    pub fn value(&self, position: Float3) -> Float {
        match self {
            Self::Gradient{rho0, origin: r0, lambda, direction: u0} => {
                rho0 * ((position - *r0).dot(*u0) / *lambda).exp()
            },
            Self::Uniform(rho) => { *rho },
        }
    }

    pub fn vertical_gradient(rho0: Float, z0: Float, lambda: Float) -> Result<Self> {
        let origin = Float3::new(0.0, 0.0, z0);
        let direction = Float3::new(0.0, 0.0, -1.0);
        Self::gradient(rho0, origin, lambda, direction)
    }
}
