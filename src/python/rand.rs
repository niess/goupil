use anyhow::{Context, Result};
use crate::numerics::{float::Float, rand::FloatRng};
use getrandom::getrandom;
use pyo3::prelude::*;
use rand::SeedableRng;
use serde_derive::{Deserialize, Serialize};
use super::numpy::{PyArray, PyScalar, ShapeArg};

#[cfg(not(feature = "f32"))]
use rand_pcg::Pcg64Mcg as Generator;
#[cfg(feature = "f32")]
use rand_pcg::Pcg32 as Generator;

#[cfg(not(feature = "f32"))]
pub(super) type Index = u128;
#[cfg(feature = "f32")]
pub(super) type Index = u64;


#[derive(Deserialize, Serialize)]
#[pyclass(name = "RandomStream", module = "goupil")]
pub struct PyRandomStream {
    pub(crate) generator: Generator,

    /// Random stream current index.
    #[pyo3(get)]
    index: Index,

    /// Random stream seed.
    #[pyo3(get)]
    seed: u128,
}

#[pymethods]
impl PyRandomStream {
    #[new]
    #[pyo3(signature = (*, seed=None, index=None))]
    pub fn new(seed: Option<u128>, index: Option<arg::Index>) -> Result<Self> {
        #[cfg(not(feature = "f32"))]
        let generator = Generator::new(0xCAFEF00DD15EA5E5);
        #[cfg(feature = "f32")]
        let generator = Generator::new(0xCAFEF00DD15EA5E5, 0xA02BDBF7BB3C0A7);
        let mut obj = Self { generator, seed: 0, index: 0 };
        obj.initialise(seed)?;
        if index.is_some() {
            obj.set_index(index)?;
        }
        Ok(obj)
    }

    #[setter]
    pub(super) fn set_index(&mut self, index: Option<arg::Index>) -> Result<()> {
        match index {
            None => self.initialise(Some(self.seed))?,
            Some(index) => {
                let index: Index = index.into();
                let delta: Index = index.wrapping_sub(self.index);
                self.generator.advance(delta);
                self.index = index;
            },
        }
        Ok(())
    }

    #[setter]
    fn set_seed(&mut self, seed: Option<u128>) -> Result<()> {
        self.initialise(seed)
    }

    /// Generates pseudo random number(s) following the Normal law.
    #[pyo3(name = "normal")]
    #[pyo3(signature = (shape=None, /))]
    fn py_normal(&mut self, py: Python, shape: Option<ShapeArg>) -> Result<PyObject> {
        self.generate(py, shape, Self::normal)
    }

    /// Generates pseudo random number(s) uniformly distributed over (0,1).
    #[pyo3(name = "uniform01")]
    #[pyo3(signature = (shape=None, /))]
    fn py_uniform01(&mut self, py: Python, shape: Option<ShapeArg>) -> Result<PyObject> {
        self.generate(py, shape, Self::uniform01)
    }
}

impl FloatRng for PyRandomStream {
    fn uniform01(&mut self) -> Float {
        self.index += 1;
        self.generator.uniform01()
    }
}

// Private interface.
impl PyRandomStream {
    #[cfg(not(feature = "f32"))]
    pub(super) fn index_2u64(&self) -> [u64; 2] {
        [
            (self.index >> 64) as u64,
            self.index as u64
        ]
    }

    #[cfg(feature = "f32")]
    pub(super) fn index(&self) -> u64 {
        self.index
    }

    fn initialise(&mut self, seed: Option<u128>) -> Result<()> {
        match seed {
            None => {
                let mut seed = [0_u8; 16];
                getrandom(&mut seed)
                    .with_context(|| "could not seed RandomEngine")?;
                self.generator = Generator::from_seed(seed);
                self.seed = u128::from_ne_bytes(seed);
            },
            Some(seed) => {
                self.seed = seed;
                let seed = u128::to_ne_bytes(seed);
                self.generator = Generator::from_seed(seed);
            },
        }
        self.index = 0;
        Ok(())
    }

    fn generate(
        &mut self,
        py: Python,
        shape: Option<ShapeArg>,
        func: fn(&mut Self) -> Float
    ) -> Result<PyObject> {
        match shape {
            None => {
                let value = func(self);
                let scalar = PyScalar::new(py, value)?;
                Ok(scalar.into_py(py))
            },
            Some(shape) => {
                let shape: Vec<usize> = shape.into();
                let n = shape.iter().product();
                let iter = (0..n).map(|_| func(self));
                let array = PyArray::<Float>::from_iter(py, &shape, iter)?;
                Ok(array.into_any().unbind())
            },
        }
    }
}

#[cfg(not(feature = "f32"))]
pub(super) mod arg {
    use pyo3::prelude::*;

    #[derive(FromPyObject)]
    pub enum Index {
        #[pyo3(transparent, annotation = "[u64;2]")]
        Array([u64; 2]),
        #[pyo3(transparent, annotation = "u128")]
        Scalar(u128),
    }

    impl From<Index> for super::Index {
        fn from(value: Index) -> Self {
            match value {
                Index::Array(value) => ((value[0] as u128) << 64) + (value[1] as u128),
                Index::Scalar(value) => value,
            }
        }
    }
}

#[cfg(feature = "f32")]
pub(super) mod arg {
    pub type Index = super::Index;
}
