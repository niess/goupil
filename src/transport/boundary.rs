use crate::numerics::float::{Float, Float3, Float3x3};
use serde_derive::{Deserialize, Serialize};


// ===============================================================================================
// External boundaries.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, Serialize)]
pub enum TransportBoundary {
    #[default]
    None,
    Box(BoxShape),
    Sector(usize),
    Sphere(SphereShape),
}

impl TransportBoundary {
    pub fn distance(&self, position: Float3, direction: Float3) -> Float {
        let distance = match self {
            Self::Box(shape) => shape.distance(position, direction),
            Self::Sphere(shape) => shape.distance(position, direction),
            _ => None,
        };
        distance.unwrap_or(Float::INFINITY)
    }

    pub fn inside(&self, position: Float3, sector: usize) -> bool {
        match self {
            Self::None => false,
            Self::Box(shape) => shape.inside(position),
            Self::Sector(index) => *index == sector,
            Self::Sphere(shape) => shape.inside(position),
        }
    }
}


// ===============================================================================================
// Generic geometry shape.
// ===============================================================================================

pub trait GeometryShape {
    /// Returns the distance to the shape, starting from *position* and along the given
    /// *direction*.
    fn distance(&self, position: Float3, direction: Float3) -> Option<Float>;

    /// Indicates wether or not the given *position* is inside the shape.
    fn inside(&self, position: Float3) -> bool;
}


// ===============================================================================================
// Box shape.
// ===============================================================================================

#[derive(Clone, Copy, Default, PartialEq, Deserialize, Serialize)]
pub struct BoxShape {
    pub center: Float3,
    pub size: Float3,
    pub rotation: Option<Float3x3>,
}

impl GeometryShape for BoxShape {
    fn distance(&self, position: Float3, direction: Float3) -> Option<Float> {
        let (v, u) = {
            let v = position - self.center;
            match self.rotation.as_ref() {
                None => (v, direction),
                Some(rotation) => (v * rotation, direction * rotation),
            }
        };

        // Intersection with an axis aligned box, avoiding (some) branchings.
        // Ref: https://tavianator.com/2011/ray_box.html
        let mut tmin = -Float::INFINITY;
        let mut tmax = Float::INFINITY;
        let mut update = |x: &Float, hx: Float, ux: &Float| {
            if *ux != 0.0 {
                let uxinv = 1.0 / ux;
                let tx1 = (-hx - x) * uxinv;
                let tx2 = (hx - x) * uxinv;
                tmin = tmin.max(tx1.min(tx2));
                tmax = tmax.min(tx1.max(tx2));
            }
        };
        update(&v.0, 0.5 * self.size.0, &u.0);
        update(&v.1, 0.5 * self.size.1, &u.1);
        update(&v.2, 0.5 * self.size.2, &u.2);
        if tmin <= 0.0 {
            if tmax > 0.0 { Some(tmax) } else { None }
        } else if tmax >= tmin {
            Some(tmin)
        } else {
            None
        }
    }

    fn inside(&self, position: Float3) -> bool {
        let r = {
            let mut r = position - self.center;
            if let Some(rotation) = self.rotation.as_ref() {
                r = r * rotation
            }
            r
        };
        (r.0.abs() < 0.5 * self.size.0) &&
        (r.1.abs() < 0.5 * self.size.1) &&
        (r.2.abs() < 0.5 * self.size.2)
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

impl GeometryShape for SphereShape {
    fn distance(&self, position: Float3, direction: Float3) -> Option<Float> {
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

    fn inside(&self, position: Float3) -> bool {
        (position - self.center).norm2() < self.radius * self.radius
    }
}


//================================================================================================
// Unit tests.
//================================================================================================
#[cfg(test)]
mod tests {
    use crate::numerics::tests::assert_float_eq;
    use super::*;

    #[test]
    fn inside_box() {
        let center = Float3::zero();
        let size = Float3::splat(2.0);
        let shape = BoxShape{ center, size, rotation: None };

        let position = Float3::splat(0.5);
        assert!(shape.inside(position));

        let position = Float3::splat(1.5);
        assert!(!shape.inside(position));
    }

    #[test]
    fn distance_box() {
        let center = Float3::zero();
        let size = Float3::splat(2.0);
        let shape = BoxShape{ center, size, rotation: None };

        let direction = Float3::new(0.0, 0.0, 1.0);
        let position = Float3::new(0.0, 0.0, -1.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, -0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.5);
        let position = Float3::new(0.0, 0.0, 0.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.0);
        let position = Float3::new(0.0, 0.0, 0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, 1.5);
        assert!(shape.distance(position, direction).is_none());

        let direction = Float3::new(0.0, 0.0, -1.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, 0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.5);
        let position = Float3::new(0.0, 0.0, 0.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.0);
        let position = Float3::new(0.0, 0.0, -0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, -1.5);
        assert!(shape.distance(position, direction).is_none());
    }

    #[test]
    fn inside_sphere() {
        let center = Float3::zero();
        let radius = 1.0;
        let shape = SphereShape{ center, radius };

        let position = Float3::splat(0.5);
        assert!(shape.inside(position));

        let position = Float3::splat(1.0);
        assert!(!shape.inside(position));
    }

    #[test]
    fn distance_sphere() {
        let center = Float3::zero();
        let radius = 1.0;
        let shape = SphereShape{ center, radius };

        let direction = Float3::new(0.0, 0.0, 1.0);
        let position = Float3::new(0.0, 0.0, -1.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, -0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.5);
        let position = Float3::new(0.0, 0.0, 0.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.0);
        let position = Float3::new(0.0, 0.0, 0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, 1.5);
        assert!(shape.distance(position, direction).is_none());

        let direction = Float3::new(0.0, 0.0, -1.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, 0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.5);
        let position = Float3::new(0.0, 0.0, 0.0);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 1.0);
        let position = Float3::new(0.0, 0.0, -0.5);
        assert_float_eq!(shape.distance(position, direction).unwrap(), 0.5);
        let position = Float3::new(0.0, 0.0, -1.5);
        assert!(shape.distance(position, direction).is_none());
    }
}
