use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};


//================================================================================================
// Floating point type (defined at compilation).
//================================================================================================
#[cfg(not(feature = "f32"))]
pub type Float = f64;
#[cfg(feature = "f32")]
pub type Float = f32;


//================================================================================================
// Vectorized floating point type.
//================================================================================================
#[derive(Copy, Clone, Default, PartialEq)]
pub struct Float3 (pub Float, pub Float, pub Float);

impl Float3 {
    // Constructors.
    pub const fn new(x: Float, y: Float, z: Float) -> Float3 { Float3(x, y, z) }
    pub const fn splat(x: Float) -> Float3 { Float3(x, x, x) }
    pub const fn zero() -> Float3 { Float3(0.0, 0.0, 0.0) }

    // L2-norm and related functions.
    pub fn norm(self) -> Float {
        self.dot(self).sqrt()
    }
    pub fn norm2(self) -> Float {
        (self.0 * self.0 + self.1 * self.1) + self.2 * self.2
    }
    pub fn normalise(&mut self) {
        *self /= self.norm();
    }
    pub fn unit(self) -> Self {
        self / self.norm()
    }

    // Inner and outer products.
    pub fn dot(self, rhs: Self) -> Float {
        (self.0 * rhs.0 + self.1 * rhs.1) + self.2 * rhs.2
    }
    pub fn cross(self, rhs: Self) -> Self {
        Self (
            self.1 * rhs.2 - self.2 * rhs.1,
            self.2 * rhs.0 - self.0 * rhs.2,
            self.0 * rhs.1 - self.1 * rhs.0,
        )
    }

    // In-place transforms.
    pub fn reverse(&mut self) {
        *self = Self (-self.0, -self.1, -self.2);
    }
    pub fn rotate(&mut self, cos_theta: Float, phi: Float) {
        // Compute (and check) the sine.
        let sin_theta = {
            let stsq = 1.0 - cos_theta * cos_theta;
            if stsq < 0.0 { return }
            stsq.sqrt()
        };

        // Get norm and unit vector.
        let norm = self.norm();
        let direction = *self / norm;

        // Generate co-vectors for the local basis.
        let e: Float3 = {
            let a0 = direction.0.abs();
            let a1 = direction.1.abs();
            let a2 = direction.2.abs();

            if a0 < a1 {
                if a0 < a2 {
                    Float3::new(1.0, 0.0, 0.0)
                } else {
                    Float3::new(0.0, 0.0, 1.0)
                }
            } else {
                if a1 < a2 {
                    Float3::new(0.0, 1.0, 0.0)
                } else {
                    Float3::new(0.0, 0.0, 1.0)
                }
            }
        };

        let mut u0 = direction.cross(e);
        u0.normalise();
        let u1 = u0.cross(direction);

        // Apply the rotation.
        *self = (norm * cos_theta) * direction +
                (norm * sin_theta) * (phi.cos() * u0 + phi.sin() * u1);
    }
}

impl fmt::Display for Float3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.0, self.1, self.2)
    }
}

macro_rules! impl_binary_operator {
    ($trait:ident, $func:ident, $op:tt) => {
        impl $trait for Float3 {
            type Output = Self;
            fn $func (self, rhs: Self) -> Self {
                Self (self.0 $op rhs.0, self.1 $op rhs.1, self.2 $op rhs.2)
            }
        }
        impl $trait<Float3> for Float {
            type Output = Float3;
            fn $func (self, rhs: Float3) -> Float3 {
                Float3 (self $op rhs.0, self $op rhs.1, self $op rhs.2)
            }
        }
        impl $trait<Float> for Float3 {
            type Output = Self;
            fn $func (self, rhs: Float) -> Self {
                Self (self.0 $op rhs, self.1 $op rhs, self.2 $op rhs)
            }
        }
    }
}

impl_binary_operator!(Add, add, +);
impl_binary_operator!(Mul, mul, *);
impl_binary_operator!(Sub, sub, -);

impl Div<Float> for Float3 {
    type Output = Self;
    fn div(self, rhs: Float) -> Self {
        let tmp = 1.0 / rhs;
        Self (self.0 * tmp, self.1 * tmp, self.2 * tmp)
    }
}

macro_rules! impl_assign_operator {
    ($trait:ident, $func:ident, $op:tt) => {
        impl $trait for Float3 {
            fn $func (&mut self, rhs: Self) {
                *self = Self (self.0 $op rhs.0, self.1 $op rhs.1, self.2 $op rhs.2);
            }
        }
        impl $trait<Float> for Float3 {
            fn $func (&mut self, rhs: Float) {
                *self = Self (self.0 $op rhs, self.1 $op rhs, self.2 $op rhs);
            }
        }
    }
}

impl_assign_operator!(AddAssign, add_assign, +);
impl_assign_operator!(MulAssign, mul_assign, *);
impl_assign_operator!(SubAssign, sub_assign, -);

impl DivAssign<Float> for Float3 {
    fn div_assign(&mut self, rhs: Float) {
        let tmp = 1.0 / rhs;
        *self = Self (self.0 * tmp, self.1 * tmp, self.2 * tmp);
    }
}

impl Neg for Float3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self (-self.0, -self.1, -self.2)
    }
}

impl From<[Float; 3]> for Float3 {
    fn from(array: [Float; 3]) -> Self {
        Self::new(array[0], array[1], array[2])
    }
}

impl From<Float3> for [Float; 3] {
    fn from(array: Float3) -> Self {
        [array.0, array.1, array.2]
    }
}

impl<'de> Deserialize<'de> for Float3 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let a: [Float; 3] = Deserialize::deserialize(deserializer)?;
        Ok(a.into())
    }
}

impl Serialize for Float3 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let a: [Float; 3] = [self.0, self.1, self.2];
        Serialize::serialize(&a, serializer)
    }
}

//================================================================================================
// Matrix floating point type.
//================================================================================================
#[derive(Copy, Clone, Default, PartialEq)]
pub struct Float3x3 ([Float; 9]);

impl From<[Float; 9]> for Float3x3 {
    fn from(a: [Float; 9]) -> Self {
        Self(a)
    }
}

impl From<[[Float; 3]; 3]> for Float3x3 {
    fn from(a: [[Float; 3]; 3]) -> Self {
        let a = [
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2],
        ];
        Self(a)
    }
}

impl From<Float3x3> for [Float; 9] {
    fn from(m: Float3x3) -> Self {
        m.0
    }
}

impl AsRef<[Float]> for Float3x3 {
    fn as_ref(&self) -> &[Float] {
        &self.0
    }
}

impl fmt::Display for Float3x3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, ..., {})", self.0[0], self.0[8])
    }
}

impl Index<(usize, usize)> for Float3x3 {
    type Output = Float;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let k = index.0 * 3 + index.1;
        &self.0[k]
    }
}

impl Mul<Float3> for &Float3x3 {
    type Output = Float3;
    fn mul (self, rhs: Float3) -> Self::Output {
        let x = self.0[0] * rhs.0 + self.0[1] * rhs.1 + self.0[2] * rhs.2;
        let y = self.0[3] * rhs.0 + self.0[4] * rhs.1 + self.0[5] * rhs.2;
        let z = self.0[6] * rhs.0 + self.0[7] * rhs.1 + self.0[8] * rhs.2;
        Float3::new(x, y, z)
    }
}

impl Mul<&Float3x3> for Float3 {
    type Output = Float3;
    fn mul (self, rhs: &Float3x3) -> Self::Output {
        let x = self.0 * rhs.0[0] + self.1 * rhs.0[3] + self.2 * rhs.0[6];
        let y = self.0 * rhs.0[1] + self.1 * rhs.0[4] + self.2 * rhs.0[7];
        let z = self.0 * rhs.0[2] + self.1 * rhs.0[5] + self.2 * rhs.0[8];
        Self::new(x, y, z)
    }
}

impl<'de> Deserialize<'de> for Float3x3 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let a: [Float; 9] = Deserialize::deserialize(deserializer)?;
        Ok(a.into())
    }
}

impl Serialize for Float3x3 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let a: &[Float] = self.as_ref();
        Serialize::serialize(&a, serializer)
    }
}


//================================================================================================
// Unit tests.
//================================================================================================
#[cfg(test)]
mod tests {
    use super::super::consts::{FRAC_1_SQRT_2, PI};
    use super::super::tests::{assert_float_eq, assert_float3_eq};
    use super::*;

    #[test]
    fn add() {
        let mut v0 = Float3::new(1.0, 2.0, 3.0);
        let mut v1 = Float3::new(4.0, 5.0, 6.0);
        const V2: Float3 = Float3::new(5.0, 7.0, 9.0);
        assert!(v0 + v1 == V2);
        assert!(v0 + 1.0 == v0 + Float3::splat(1.0));
        assert!(1.0 + v0 == Float3::splat(1.0) + v0);

        v0 += v1;
        assert!(v0 == V2);
        v1 += 1.0;
        assert!(v1 == Float3::new(5.0, 6.0, 7.0));
    }

    #[test]
    fn constexpr() {
        const V0: Float3 = Float3::new(1.0, 2.0, 3.0);
        let mut v1 = Float3::new(V0.0, V0.1, V0.2);
        assert!(v1 == V0);

        const V2: Float3 = Float3::splat(1.0);
        const V3: Float3 = Float3::zero();
        v1 += 1.0;
        assert!(v1 == V0 + V2 + V3);
    }

    #[test]
    fn copy() {
        let v0 = Float3::new(1.0, 2.0, 3.0);
        let v1 = v0;
        assert!(v0 == v1);
    }

    #[test]
    fn div() {
        let mut v0 = Float3::new(2.0, 4.0, 6.0);
        const V1: Float3 = Float3::new(1.0, 2.0, 3.0);
        assert!(v0 / 2.0 == V1);
        v0 /= 2.0;
        assert!(v0 == V1);
    }

    #[test]
    fn equal() {
        let v0 = Float3::new(1.0, 2.0, 3.0);
        let v1 = Float3::new(1.0, 2.0, 3.0);
        assert!(v0 == v1);
        assert!(v0 != Float3::zero());
    }

    #[test]
    fn mul() {
        let mut v0 = Float3::new(1.0, 2.0, 3.0);
        let mut v1 = Float3::new(4.0, 5.0, 6.0);
        const V2: Float3 = Float3::new(4.0, 10.0, 18.0);
        assert!(v0 * v1 == V2);
        assert!(v0 * 2.0 == v0 * Float3::splat(2.0));
        assert!(2.0 * v0 == Float3::splat(2.0) * v0);

        v0 *= v1;
        assert!(v0 == V2);
        v1 *= 2.0;
        assert!(v1 == Float3::new(8.0, 10.0, 12.0));
    }

    #[test]
    fn neg() {
        let v0 = Float3::new(1.0, 2.0, 3.0);
        assert!(-v0 == v0 * -1.0);
    }

    #[test]
    fn norm() {
        const V: Float3 = Float3::new(1.0, 2.0, 3.0);
        const F: Float = 14.0;
        assert_eq!(V.norm2(), F);
        assert_eq!(V.norm(), F.sqrt());

        let u = V.unit();
        assert_float_eq!(u.norm(), 1.0);

        let mut v = V;
        v.normalise();
        assert_float_eq!(v.norm(), 1.0);
    }

    #[test]
    fn product() {
        let v0 = Float3::new(1.0, 1.0, 1.0);
        let v1 = Float3::new(0.0, 1.0, 1.0);
        let v2 = v0.cross(v1);
        assert_eq!(v2.dot(v0), 0.0);
        assert_eq!(v2.dot(v1), 0.0);
        let cos_theta = v0.dot(v1) / (v0.norm() * v1.norm());
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        assert_float_eq!(v2.norm(), v0.norm() * v1.norm() * sin_theta);
    }

    #[test]
    fn rotate() {
        let mut v0 = Float3::new(1.0, 0.0, 1.0);
        let mut v1 = v0;
        v0.rotate(-1.0, 0.0);
        v1.reverse();
        assert_float3_eq!(v0, v1);

        v1.reverse();
        v1.rotate(0.0, PI / 2.0);
        assert_eq!(v0.dot(v1), 0.0);

        let v2 = v1;
        v1.rotate(0.5, PI / 4.0);
        assert_float_eq!(v2.dot(v1), 1.0);

        let mut v3 = Float3::new(0.0, 0.0, 1.0);
        v3.rotate(0.0, PI / 4.0);
        const TMP: Float = FRAC_1_SQRT_2;
        assert_float_eq!(v3.0, -TMP);
        assert_float_eq!(v3.1, TMP);
    }

    #[test]
    fn sub() {
        let mut v0 = Float3::new(1.0, 2.0, 3.0);
        let mut v1 = Float3::new(4.0, 6.0, 8.0);
        const V2: Float3 = Float3::new(-3.0, -4.0, -5.0);
        assert!(v0 - v1 == V2);
        assert!(v0 - 1.0 == v0 - Float3::splat(1.0));
        assert!(1.0 - v0 == Float3::splat(1.0) - v0);

        v0 -= v1;
        assert!(v0 == V2);
        v1 -= 1.0;
        assert!(v1 == Float3::new(3.0, 5.0, 7.0));
    }
}
