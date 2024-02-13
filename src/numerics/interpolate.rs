use serde_derive::{Deserialize, Serialize};
use std::ops::{Index, IndexMut};
use super::float::Float;
use super::grids::{Grid, GridCoordinate};


//================================================================================================
// Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).
//
// The derivative coefficients are computed using the method of Fritsch and Butland. For boundary
// conditions a 3 points finite difference is used.
//
// Note: The implementation assumes that inputs always have consistent sizes. It will panic
// otherwise. Note also that function values at grid nodes are not managed by the interpolator
// object. Only, derivative coefficients are handled.
//
// References:
//  F. N. Fristch and J. Butland, SIAM J. Sci. Stat. Comput. (1984)
//
//================================================================================================

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct CubicInterpolator (Vec<Float>);

impl CubicInterpolator {
    pub fn new(n: usize) -> Self {
        Self (vec![0.0; n])
    }

    // Initialise derivative coefficients.
    pub fn initialise<T: Grid>(&mut self, x: &T, y: &[Float], der: bool) {
        let n = x.len();
        if n == 0 {
            return;
        } else if n == 1 {
            if !der { self.0[0] = 0.0 }
        } else if n == 2 {
            if !der {
                let d = (y[1] - y[0]) / x.width(0);
                self.0[0] = d;
                self.0[1] = d;
            }
        } else {
            for i in 1..(n - 1) {
                let h1 = x.width(i - 1);
                if h1 == 0.0 {
                    self.0[i] = 0.0;
                    continue;
                }
                if der {
                    let a = self.0[i] / h1;
                    let b = self.0[i + 1] / h1;
                    let c = 2.0 * a + b - 3.0;
                    if 3.0 * a * (a + b - 2.0) >= c * c { continue }
                }

                let h2 = x.width(i);
                if h2 == 0.0 {
                    self.0[i] = 0.0;
                    continue;
                }
                let s1 = (y[i] - y[i - 1]) / h1;
                let s2 = (y[i + 1] - y[i]) / h2;

                let tmp = s1 * s2;
                if tmp > 0.0 {
                    let a = (h1 + 2.0 * h2) / (3.0 * (h1 + h2));
                    self.0[i] = tmp / ((1.0 - a) * s1 + a * s2);
                } else {
                    self.0[i] = 0.0;
                }
            }

            if !der {
                let diff3 = |h1: Float, h2: Float, y0, y1, y2| -> Float {
                    let h2 = h1 + h2;
                    let delta = h1 * h2 * (h2 - h1);
                    if delta == 0.0 { return 0.0 }
                    let c1 = h2 * h2 / delta;
                    let c2 = -h1 * h1 / delta;
                    let c0 = -(c1 + c2);
                    c0 * y0 + c1 * y1 + c2 * y2
                };

                self.0[0] = diff3(x.width(0), x.width(1), y[0], y[1], y[2]);
                self.0[n - 1] = diff3(
                    -x.width(n - 2), -x.width(n - 3), y[n - 1], y[n - 2], y[n - 3]);
            }
        }
    }

    pub fn interpolate<T: Grid>(&self, grid: &T, y: &[Float], x: Float) -> Option<Float> {
        match grid.transform(x) {
            GridCoordinate::Inside(i, d) => {
                Some(self.interpolate_raw(i, d, grid.width(i), y))
            },
            _ => None,
        }
    }

    pub fn interpolate_raw(&self, i: usize, d: Float, dx: Float, y: &[Float]) -> Float {
        let p0 = y[i];
        let p1 = y[i + 1];
        let m0 = self.0[i] * dx;
        let m1 = self.0[i + 1] * dx;
        let c2 = -3.0 * (p0 - p1) - 2.0 * m0 - m1;
        let c3 = 2.0 * (p0 - p1) + m0 + m1;
        p0 + d * (m0 + d * (c2 + d * c3))
    }
}

impl Index<usize> for CubicInterpolator {
    type Output = Float;

    fn index(&self, i: usize) -> &Self::Output { &self.0[i] }
}

impl IndexMut<usize> for CubicInterpolator {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output { &mut self.0[i] }
}

impl AsMut<Vec<Float>> for CubicInterpolator {
    fn as_mut(&mut self) -> &mut Vec<Float> { &mut self.0 }
}

impl AsRef<Vec<Float>> for CubicInterpolator {
    fn as_ref(&self) -> &Vec<Float> { &self.0 }
}


//================================================================================================
// A bi-linear interpolator.
//================================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct BilinearInterpolator {
    columns: usize,
    table: Vec<Float>,
}

impl BilinearInterpolator {
    pub fn interpolate<GI: Grid, GJ: Grid>(&self, grid_x: &GJ, x: Float, grid_y: &GI, y: Float)
        -> Float {

        let (i, hi) = grid_y.transform(y).clamp();
        let (j, hj) = grid_x.transform(x).clamp();
        self.interpolate_raw(i, hi, j, hj)
    }

    pub fn interpolate_raw(&self, i: usize, hi: Float, j: usize, hj: Float) -> Float {
        let f00 = self[(i, j)];
        if f00.is_nan() { return f00 }
        let f01 = self[(i, j + 1)];
        if f01.is_nan() { return f01 }
        let f10 = self[(i + 1, j)];
        if f10.is_nan() { return f10 }
        let f11 = self[(i + 1, j + 1)];
        if f11.is_nan() { return f11 }
        let ti = 1.0 - hi;
        let tj = 1.0 - hj;
        f00 * ti * tj + f01 * ti * hj + f10 * hi * tj + f11 * hi * hj
    }

    pub fn new(lines: usize, columns: usize) -> Self {
        let table = vec![0.0; lines * columns];
        Self { columns, table }
    }

    pub fn shape(&self) -> (usize, usize) {
        let lines = self.table.len() / self.columns;
        (lines, self.columns)
    }
}

impl Index<(usize, usize)> for BilinearInterpolator {
    type Output = Float;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let k = index.0 * self.columns + index.1;
        &self.table[k]
    }
}

impl IndexMut<(usize, usize)> for BilinearInterpolator {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let k = index.0 * self.columns + index.1;
        &mut self.table[k]
    }
}

impl AsRef<Vec<Float>> for BilinearInterpolator {
    fn as_ref(&self) -> &Vec<Float> { &self.table }
}


//================================================================================================
// Unit tests.
//================================================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use super::super::grids::LogGrid;
    use super::super::tests::assert_float_eq;

    #[test]
    fn cubic_interpolator() {
        let xg = LogGrid::new(1.0, 10.0, 901);
        let f = |x: Float| -> Float { x * x * x - 2.0 * x * x + 1.0 };
        let yg: Vec<Float> = (0..xg.len()).map(|j| {f(xg[j])}).collect();
        let mut i = CubicInterpolator::new(xg.len());
        i.initialise(&xg, &yg, false);
        let eps = match std::mem::size_of::<Float>() {
            4 => 1e-4,
            _ => 1e-5,
        };
        for j in 0..901 {
            let x: Float = (j as Float) * 0.01 + 1.0;
            let y = f(x);
            let dy = eps * y.abs();
            assert_float_eq!(y, i.interpolate(&xg, &yg, x).unwrap_or(0.0), dy);
        }
    }
}
