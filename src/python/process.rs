use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::physics::process::absorption::AbsorptionMode;
use crate::physics::process::compton::{
    self,
    ComptonModel,
    compute::ComptonComputer,
    sample::ComptonSampler,
    ComptonMode,
    ComptonMethod,
};
use crate::physics::process::rayleigh::{RayleighMode, sample::RayleighSampler};
use crate::transport::TransportSettings;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use super::macros::{key_error, not_implemented_error, value_error};
use super::materials::MaterialLike;
use super::numpy::{ArrayOrFloat, PyArray, PyArrayMethods, PyScalar};
use super::rand::PyRandomStream;


// ===============================================================================================
// Python wrapper for Compton process.
// ===============================================================================================
#[pyclass(name = "ComptonProcess", module = "goupil")]
pub struct PyComptonProcess {
    computer: ComptonComputer,
    sampler: ComptonSampler,
}

#[pymethods]
impl PyComptonProcess {

    #[getter]
    fn get_method(&self) -> &str {
        self.sampler.method.into()
    }

    #[setter]
    fn set_method(&mut self, value: &str) -> Result<()> {
        self.sampler.method = ComptonMethod::try_from(value)?;
        Ok(())
    }

    #[getter]
    fn get_mode(&self) -> Option<&str> {
        match self.computer.mode {
            ComptonMode::None => None,
            _ => Some(self.computer.mode.into())
        }
    }

    #[setter]
    fn set_mode(&mut self, value: Option<&str>) -> Result<()> {
        match value {
            None => {
                self.computer.mode = ComptonMode::None;
                self.sampler.mode = ComptonMode::None;
            },
            Some(value) => {
                let value = ComptonMode::try_from(value)?;
                self.computer.mode = value;
                self.sampler.mode = value;
            }
        }
        Ok(())
    }

    #[getter]
    fn get_model(&self) -> &str {
        self.computer.model.into()
    }

    #[setter]
    fn set_model(&mut self, value: &str) -> Result<()> {
        let value = ComptonModel::try_from(value)?;
        self.computer.model = value;
        self.sampler.model = value;
        Ok(())
    }

    #[getter]
    fn get_precision(&self) -> Float {
        self.computer.precision
    }

    #[setter]
    fn set_precision(&mut self, value: Float) -> Result<()> {
        if value <= 0.0 {
            value_error!(
                "bad precision (expected a positive value, found {})",
                value
            )
        }
        self.computer.precision = value;
        Ok(())
    }

    #[new]
    #[pyo3(signature = (**kwargs))]
    fn new(kwargs: Option<&Bound<PyDict>>) -> Result<Self> {
        let mut method = ComptonMethod::default();
        let mut mode = ComptonMode::default();
        let mut model = ComptonModel::default();
        let mut precision: Option<Float> = None;
        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key: String = key.extract()?;
                match key.as_str() {
                    "method" => {
                        let value: String = value.extract()?;
                        method = ComptonMethod::try_from(value.as_str())?;
                    },
                    "mode" => {
                        let value: String = value.extract()?;
                        mode = ComptonMode::try_from(value.as_str())?;
                    },
                    "model" => {
                        let value: String = value.extract()?;
                        model = ComptonModel::try_from(value.as_str())?;
                    },
                    "precision" => {
                        let value: Float = value.extract()?;
                        precision = Some(value);
                    },
                    _ => key_error!(
                        "bad attribute (expected one of 'method', 'mode', 'model', 'precision', 
                            found '{}'",
                        key
                    ),
                }
            }
        }
        if let Err(err) = compton::validate(model, mode, method) {
            not_implemented_error!("{}", err)
        }
        let computer = ComptonComputer::new(model, mode);
        let sampler = ComptonSampler::new(model, mode, method);
        let mut object = Self { computer, sampler };
        if let Some(precision) = precision {
            object.set_precision(precision)?;
        }
        Ok(object)
    }

    fn __repr__(&self) -> String {
        let mut s = String::from("ComptonProcess(");
        let prefixes = vec!["", ", "];
        let mut prefix_index = 0;
        if self.sampler.method != ComptonMethod::default() {
            let value: &str = self.sampler.method.into();
            s.push_str(&format!(
                "method=\"{}\"",
                value
            ));
            prefix_index += 1;
        }
        if self.sampler.mode != ComptonMode::default() {
            let value: &str = self.sampler.mode.into();
            s.push_str(&format!(
                "{}mode=\"{}\"",
                prefixes[prefix_index],
                value
            ));
            if prefix_index == 0 {
                prefix_index = 1;
            }
        }
        if self.sampler.model != ComptonModel::default() {
            let value: &str = self.sampler.model.into();
            s.push_str(&format!(
                "{}model=\"{}\"",
                prefixes[prefix_index],
                value
            ));
            if prefix_index == 0 {
                prefix_index = 1;
            }
        }
        if self.computer.precision != 1.0 {
            s.push_str(&format!(
                "{}precision={}",
                prefixes[prefix_index],
                self.computer.precision
            ));
        }
        s.push_str(")");
        s
    }

    #[pyo3(signature = (energy, material, *, energy_min=None, energy_max=None))]
    fn cross_section(
        &self,
        py: Python,
        energy: ArrayOrFloat,
        material: MaterialLike,
        energy_min: Option<Float>,
        energy_max: Option<Float>
    ) -> Result<PyObject> {
        let electrons = material.get_electrons()?;
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let result = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = self.computer.cross_section(
                        energy.get(i)?,
                        energy_min,
                        energy_max,
                        &electrons,
                    )?;
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy) => {
                let result = self.computer.cross_section(
                    energy,
                    energy_min,
                    energy_max,
                    &electrons,
                )?;
                let result = PyScalar::<Float>::new(py, result)?;
                result.into_py(py)
            }
        };
        Ok(result)
    }

    fn dcs(
        &self,
        py: Python,
        energy_in: Float,
        energy_out: ArrayOrFloat,
        material: MaterialLike
    ) -> Result<PyObject> {
        if let ComptonModel::Penelope = self.computer.model {
            not_implemented_error!(
                "DCS w.r.t. energy is not implemented for Penelope model"
            );
        }
        let electrons = material.get_electrons()?;
        let result: PyObject = match energy_out {
            ArrayOrFloat::Array(energy_out) => {
                let result = PyArray::<Float>::empty(py, &energy_out.shape())?;
                let n = energy_out.size();
                for i in 0..n {
                    let v = self.computer.dcs(
                        energy_in,
                        energy_out.get(i)?,
                        &electrons,
                    );
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy_out) => {
                let result = self.computer.dcs(
                    energy_in,
                    energy_out,
                    &electrons,
                );
                result.into_py(py)
            }
        };
        Ok(result)
    }

    fn dcs_support(&self, py: Python, energy: ArrayOrFloat) -> Result<PyObject> {
        if let ComptonModel::Penelope = self.computer.model {
            not_implemented_error!(
                "DCS support w.r.t. energy is not implemented for Penelope model"
            );
        }
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let energy_min = PyArray::<Float>::empty(py, &energy.shape())?;
                let energy_max = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let (min, max) = self.computer.dcs_support(energy.get(i)?);
                    energy_min.set(i, min)?;
                    energy_max.set(i, max)?;
                }
                let energy_min: PyObject = energy_min.into_py(py);
                let energy_max: PyObject = energy_max.into_py(py);
                (energy_min, energy_max).into_py(py)
            },
            ArrayOrFloat::Float(energy) => {
                let result = self.computer.dcs_support(energy);
                result.into_py(py)
            }
        };
        Ok(result)
    }

    #[pyo3(signature = (energy, material, *, rng=None))]
    fn sample<'py>(
        &self,
        py: Python<'py>,
        energy: ArrayOrFloat,
        material: MaterialLike,
        rng: Option<&Bound<'py, PyRandomStream>>,
    )
    -> Result<PyObject> {
        // Prepare material and generator.
        // let material = material.get(py)?;
        let mut registry;
        let material = match &material {
            MaterialLike::Record(record) => record.get(py)?,
            _ => {
                let definition = material.unpack()?;
                registry = MaterialLike::new_registry(py, &definition)?;
                let mut settings = TransportSettings::default();
                settings.compton_method = self.sampler.method;
                settings.compton_mode = self.sampler.mode;
                settings.compton_model = self.sampler.model;
                settings.absorption = AbsorptionMode::None;
                settings.rayleigh = RayleighMode::None;
                registry.compute(&settings, None, None, None)?;
                registry.get(definition.name())?
            }
        };

        let rng = match rng {
            None => Bound::new(py, PyRandomStream::new(None, None)?)?,
            Some(rng) => rng.clone(),
        };
        let mut rng = rng.borrow_mut();

        // Generate samples.
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let energy_out = PyArray::<Float>::empty(py, &energy.shape())?;
                let cos_theta = PyArray::<Float>::empty(py, &energy.shape())?;
                let weight = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let momentum_in = Float3::new(0.0, 0.0, energy.get(i)?);
                    let sample = self.sampler.sample(
                        &mut rng.generator,
                        momentum_in,
                        material,
                        None,
                    )?;
                    let e = sample.momentum_out.norm();
                    energy_out.set(i, e)?;
                    cos_theta.set(i, sample.momentum_out.2 / e)?;
                    weight.set(i, sample.weight)?;
                }
                let energy_out: PyObject = energy_out.into_py(py);
                let cos_theta: PyObject = cos_theta.into_py(py);
                let weight: PyObject = weight.into_py(py);
                (energy_out, cos_theta, weight).into_py(py)
            },
            ArrayOrFloat::Float(energy) => {
                let momentum_in = Float3::new(0.0, 0.0, energy);
                let sample = self.sampler.sample(
                    &mut rng.generator,
                    momentum_in,
                    material,
                    None,
                )?;
                let energy_out = sample.momentum_out.norm();
                let cos_theta = sample.momentum_out.2 / energy_out;
                (energy_out, cos_theta, sample.weight).into_py(py)
            }
        };
        Ok(result)
    }
}


// ===============================================================================================
// Python wrapper for Rayleigh process.
// ===============================================================================================
#[pyclass(name = "RayleighProcess", module = "goupil")]
pub struct PyRayleighProcess ();

#[pymethods]
impl PyRayleighProcess {
    #[staticmethod]
    fn cross_section(
        py: Python,
        energy: ArrayOrFloat,
        material: MaterialLike,
    ) -> Result<PyObject> {
        let cross_section = material.rayleigh_cross_section(py)?;
        let compute = |energy: Float| -> Float {
            cross_section.interpolate(energy)
        };
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let result = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = compute(energy.get(i)?);
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy) => compute(energy).into_py(py),
        };
        Ok(result)
    }

    #[staticmethod]
    fn dcs(
        py: Python,
        energy: Float,
        cos_theta: ArrayOrFloat,
        material: MaterialLike,
    ) -> Result<PyObject> {
        let sampler = RayleighSampler::new(RayleighMode::FormFactor);
        let form_factor = material.rayleigh_form_factor(py)?;
        let result: PyObject = match cos_theta {
            ArrayOrFloat::Array(cos_theta) => {
                let result = PyArray::<Float>::empty(py, &cos_theta.shape())?;
                let n = cos_theta.size();
                for i in 0..n {
                    let v = sampler.dcs(energy, cos_theta.get(i)?, &form_factor)?;
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(cos_theta) => {
                let result = sampler.dcs(energy, cos_theta, &form_factor)?;
                result.into_py(py)
            },
        };
        Ok(result)
    }

    #[staticmethod]
    #[pyo3(signature = (energy, material, *, rng=None))]
    fn sample<'py>(
        py: Python<'py>,
        energy: ArrayOrFloat,
        material: MaterialLike,
        rng: Option<&Bound<'py, PyRandomStream>>,
    )
    -> Result<PyObject> {
        let sampler = RayleighSampler::new(RayleighMode::FormFactor);
        let rng = match rng {
            None => Bound::new(py, PyRandomStream::new(None, None)?)?,
            Some(rng) => rng.clone(),
        };
        let mut rng = rng.borrow_mut();
        let form_factor = material.rayleigh_form_factor(py)?;
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let result = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = sampler.sample_angle(
                        &mut rng.generator, energy.get(i)?, &form_factor)?;
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy) => {
                let result = sampler.sample_angle(&mut rng.generator, energy, &form_factor)?;
                result.into_py(py)
            },
        };
        Ok(result)
    }
}


// ===============================================================================================
// Python wrapper for absorption process.
// ===============================================================================================
#[pyclass(name = "AbsorptionProcess", module = "goupil")]
pub struct PyAbsorptionProcess ();

#[pymethods]
impl PyAbsorptionProcess {
    #[staticmethod]
    fn cross_section(
        py: Python,
        energy: ArrayOrFloat,
        material: MaterialLike,
    ) -> Result<PyObject> {
        let cross_section = material.absorption_cross_section(py)?;
        let compute = |energy: Float| -> Float {
            cross_section.interpolate(energy)
        };
        let result: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let result = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = compute(energy.get(i)?);
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy) => compute(energy).into_py(py),
        };
        Ok(result)
    }
}
