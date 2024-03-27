use serde_derive::{Deserialize, Serialize};
use std::ops::{Index, IndexMut};
use super::float::Float;


//================================================================================================
// Interface for one-dimensional grids.
//================================================================================================

pub enum GridCoordinate {
    Above(usize),
    Below,
    Inside(usize, Float),
}

impl GridCoordinate {
    pub fn clamp(&self) -> (usize, Float) {
        match self {
            GridCoordinate::Above(n) => (n - 2, 1.0),
            GridCoordinate::Below => (0, 0.0),
            GridCoordinate::Inside(i, d) => (*i, *d),
        }
    }
}

pub trait Grid {
    // Grid length, i.e. the number of nodes.
    fn len(&self) -> usize;

    // Transform to grid coordinates.
    fn transform(&self, x: Float) -> GridCoordinate;

    // Width of the i^th grid segment.
    fn width(&self, i: usize) -> Float;
}


//================================================================================================
// Linear grid, i.e. with a constant spacing.
//================================================================================================

#[derive(Clone, Default, Deserialize, PartialEq, Serialize)]
pub struct LinearGrid {
    xmin: Float,
    xmax: Float,
    dx: Float,
    n: usize,
    i: usize,
}

impl LinearGrid {
    pub fn get(&self, i: usize) -> Float {
        if i == 0 {
            self.xmin
        } else if i == self.n - 1 {
            self.xmax
        } else {
            self.xmin + (i as Float) * self.dx
        }
    }

    pub fn iter(&self) -> Self {
        let mut clone = self.clone();
        clone.i = 0;
        clone
    }

    pub fn new(xmin: Float, xmax: Float, n: usize) -> Self {
        let dx = (xmax - xmin) / (n - 1) as Float;
        let i: usize = 0;
        Self { xmin, xmax, dx, n, i }
    }
}

impl Into<Vec<Float>> for LinearGrid {
    fn into(self) -> Vec<Float> {
        (0..self.n).map(|i: usize| self.get(i)).collect()
    }
}

impl Grid for LinearGrid {
    fn len(&self) -> usize { self.n }

    fn transform(&self, x: Float) -> GridCoordinate {
        let tmp = (x - self.xmin) / self.dx;
        if tmp < 0.0 {
            return GridCoordinate::Below
        }
        let i = tmp as usize;
        if i > self.n - 1 {
            GridCoordinate::Above(self.n)
        } else if i == self.n - 1 {
            if x <= self.xmax {
                GridCoordinate::Inside(self.n - 2, 1.0)
            } else {
                GridCoordinate::Above(self.n)
            }
        } else {
            let h = tmp - (i as Float);
            GridCoordinate::Inside(i, h)
        }
    }

    fn width(&self, _: usize) -> Float {
        self.dx
    }
}

impl Iterator for LinearGrid {
    type Item = Float;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.n {
            let x = self.get(self.i);
            self.i += 1;
            Some(x)
        } else {
            None
        }
    }
}


//================================================================================================
// Log-regular grid, i.e. with a logarithmic spacing.
//================================================================================================

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct LogGrid {
    x: Vec<Float>,
    dx: Float,
}

impl LogGrid {
    pub fn as_ptr(&self) -> *const Float {
        self.x.as_ptr()
    }

    pub fn new(xmin: Float, xmax: Float, n: usize) -> Self {
        let dx = (xmax / xmin).ln() / (n - 1) as Float;
        let mut x = vec![0.0; n];
        x[0] = xmin;
        for i in 1..(n - 1) {
            x[i] = xmin * ((i as Float) * dx).exp();
        }
        x[n - 1] = xmax;
        Self { x, dx }
    }
}

impl PartialEq for LogGrid {
    fn eq(&self, other: &Self) -> bool {
        if self.dx != other.dx {
            return false
        } else {
            let n = self.x.len();
            if other.x.len() != n {
                return false
            } else if n == 0 {
                return true
            } else {
                return self.x[0] == other.x[0] && self.x[n - 1] == other.x[n - 1]
            }
        }
    }
}

impl Index<usize> for LogGrid {
    type Output = Float;

    fn index(&self, i: usize) -> &Self::Output { &self.x[i] }
}

impl AsRef<Vec<Float>> for LogGrid {
    fn as_ref(&self) -> &Vec<Float> { &self.x }
}

impl From<LogGrid> for Vec<Float> {
    fn from(g: LogGrid) -> Self { g.x }
}

impl Grid for LogGrid {
    fn len(&self) -> usize { self.x.len() }

    fn transform(&self, x: Float) -> GridCoordinate {
        let tmp = (x / self.x[0]).ln() / self.dx;
        if tmp < 0.0 {
            return GridCoordinate::Below;
        }
        let n = self.x.len();
        let i = tmp as usize;
        if i > n - 1 {
            return GridCoordinate::Above(n);
        } else if i == n - 1 {
            if x <= self.x[n - 1] {
                return GridCoordinate::Inside(n - 2, 1.0);
            } else {
                return GridCoordinate::Above(n);
            }
        } else {
            let h = (x - self.x[i]) / (self.x[i + 1] - self.x[i]);
            return GridCoordinate::Inside(i, h);
        }
    }

    fn width(&self, i: usize) -> Float { self.x[i + 1] - self.x[i] }
}


//================================================================================================
// Unstructured grid.
//================================================================================================

#[derive(Clone, Default, Deserialize, PartialEq, Serialize)]
pub struct UnstructuredGrid(Vec<Float>);

impl UnstructuredGrid {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn new(n: usize) -> Self {
        Self ( vec![0.0; n] )
    }
}

impl Index<usize> for UnstructuredGrid {
    type Output = Float;

    fn index(&self, i: usize) -> &Self::Output { &self.0[i] }
}

impl IndexMut<usize> for UnstructuredGrid {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output { &mut self.0[i] }
}

impl AsMut<Vec<Float>> for UnstructuredGrid {
    fn as_mut(&mut self) -> &mut Vec<Float> { &mut self.0 }
}

impl AsRef<Vec<Float>> for UnstructuredGrid {
    fn as_ref(&self) -> &Vec<Float> { &self.0 }
}

impl From<Vec<Float>> for UnstructuredGrid {
    fn from(v: Vec<Float>) -> Self { Self (v) }
}

impl From<UnstructuredGrid> for Vec<Float> {
    fn from(g: UnstructuredGrid) -> Self { g.0 }
}

impl<const N: usize> From<[Float; N]> for UnstructuredGrid {
    fn from(v: [Float; N]) -> Self { Self (v.to_vec()) }
}

impl Grid for UnstructuredGrid {
    fn len(&self) -> usize { self.0.len() }

    fn transform(&self, x: Float) -> GridCoordinate {
        let mut i0 = 0;
        if x < self.0[i0] { return GridCoordinate::Below }
        let mut i1 = self.0.len() - 1;
        if x > self.0[i1] { return GridCoordinate::Above(self.0.len()) }
        else if x == self.0[i1] { return GridCoordinate::Inside(i1 - 1, 1.0) }
        while i1 > i0 + 1 {
            let i2 = (i0 + i1) / 2;
            if x < self.0[i2] {
                i1 = i2;
            } else {
                i0 = i2;
            }
        }
        let x0 = self.0[i0];
        let t = (x - x0) / (self.0[i1] - x0);
        GridCoordinate::Inside(i0, t)
    }

    fn width(&self, i: usize) -> Float { self.0[i + 1] - self.0[i] }
}


//================================================================================================
// Unit tests.
//================================================================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unstructured_grid() {
        let grid = UnstructuredGrid::from([1.0, 3.0, 7.0, 15.0]);
        assert_eq!(grid[0], 1.0);
        assert_eq!(grid.len(), 4);
        assert_eq!(grid.width(0), 2.0);
        assert_eq!(grid.transform(1.0).clamp(), (0, 0.0));
        assert_eq!(grid.transform(16.0).clamp(), (2, 1.0));
        assert_eq!(grid.transform(0.0).clamp(), (0, 0.0));
        assert_eq!(grid.transform(15.0).clamp(), (2, 1.0));
        assert_eq!(grid.transform(5.0).clamp(), (1, 0.5));
    }
}
