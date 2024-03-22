use crate::numerics::float::{Float, Float3};
use serde_derive::{Deserialize, Serialize};


// ===============================================================================================
// External boundaries.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, Serialize)]
pub enum TransportBoundary {
    #[default]
    None,
    Sector(usize),
    Sphere(SphereShape),
}

impl TransportBoundary {
    pub fn distance(&self, position: Float3, direction: Float3) -> Float {
        let distance = match self {
            Self::Sphere(sphere) => sphere.distance(position, direction),
            _ => None,
        };
        distance.unwrap_or(Float::INFINITY)
    }

    pub fn inside(&self, position: Float3, sector: usize) -> bool {
        match self {
            Self::None => false,
            Self::Sector(index) => *index == sector,
            Self::Sphere(sphere) => sphere.inside(position),
        }
    }
}


// ===============================================================================================
// Spherical shape.
// ===============================================================================================

#[derive(Clone, Copy, Default, PartialEq, Deserialize, Serialize)]
pub struct SphereShape {
    pub center: Float3,
    pub radius: Float,
}

impl SphereShape {
    pub fn distance(&self, position: Float3, direction: Float3) -> Option<Float> {
        let v = self.center - position;
        let vu = v.dot(direction);
        let h2 = v.norm2() - vu * vu;
        let r2 = self.radius * self.radius;
        if h2 > r2 {
            // No intersection case.
            None
        } else if h2 == r2 {
            // Tangent intersection case.
            if vu > 0.0 {
                Some(vu)
            } else {
                None
            }
        } else {
            // Two intersections case.
            let delta = (r2 - h2).sqrt();
            let d0 = vu + delta;
            if d0 > 0.0 {
                let d1 = vu - delta;
                if d1 > 0.0 {
                    Some(d1)
                } else {
                    Some(d0)
                }
            } else {
                None
            }
        }
    }

    pub fn inside(&self, position: Float3) -> bool {
        (position - self.center).norm2() < self.radius * self.radius
    }
}
