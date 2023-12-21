use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::physics::materials::electronic::ElectronicStructure;
use crate::physics::process::compton::{
    ComptonModel,
    compute::ComptonComputer,
    sample::ComptonSampler,
    ComptonMode,
    ComptonMethod,
};
use crate::physics::process::rayleigh::{RayleighMode, sample::RayleighSampler};
use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use pyo3::types::PyDict;
use super::macros::{key_error, value_error};
use super::materials::PyMaterialRecord;


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
    fn new(kwargs: Option<&PyDict>) -> Result<Self> {
        let mut method = ComptonMethod::default();
        let mut mode = ComptonMode::default();
        let mut model = ComptonModel::default();
        let mut precision: Option<Float> = None;
        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key: &str = key.extract()?;
                match key {
                    "method" => {
                        let value: &str = value.extract()?;
                        method = ComptonMethod::try_from(value)?;
                    },
                    "mode" => {
                        let value: &str = value.extract()?;
                        mode = ComptonMode::try_from(value)?;
                    },
                    "model" => {
                        let value: &str = value.extract()?;
                        model = ComptonModel::try_from(value)?;
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

    fn cross_section(
        &self,
        energy: Float,
        material: PyRef<PyMaterialRecord>,
        energy_min: Option<Float>,
        energy_max: Option<Float>
    ) -> Result<Float> {
        let electrons = Self::get_electrons(material.py(), &material)?;
        self.computer.cross_section(
            energy,
            energy_min,
            energy_max,
            electrons,
        )
    }

    fn dcs(
        &self,
        energy_in: Float,
        energy_out: Float,
        material: PyRef<PyMaterialRecord>
    ) -> Result<Float> {
        let electrons = Self::get_electrons(material.py(), &material)?;
        Ok(self.computer.dcs(
            energy_in,
            energy_out,
            electrons,
        ))
    }

    fn dcs_support(&self, energy_in: Float) -> (Float, Float) {
        self.computer.dcs_support(energy_in)
    }

    fn sample(
        &self,
        energy_in: Float,
        material: PyRef<PyMaterialRecord>
    )
    -> Result<(Float, Float, Float)> {
        // XXX Use PyRandomStream?
        // XXX Vectorize this method?

        // Get / format inputs.
        let mut rng = rand::thread_rng();
        let momentum_in = Float3::new(0.0, 0.0, energy_in);

        // Generate a sample.
        let py = material.py();
        let sample = self.sampler.sample(
            &mut rng,
            momentum_in,
            material.get(py)?,
            None,
        )?;

        // Format outputs.
        let energy_out = sample.momentum_out.norm();
        let cos_theta = sample.momentum_out.2 / energy_out;
        Ok((energy_out, cos_theta, sample.weight))
    }
}

// Private interface.
impl PyComptonProcess {
    fn get_electrons<'py>(
        py: Python<'py>,
        material: &'py PyMaterialRecord
    ) -> Result<&'py ElectronicStructure> {
        Ok(material
            .get(py)?
            .electrons()
            .ok_or_else(|| PyTypeError::new_err(
                "missing electronic structure (expected Some(ElectronicStructure), found None)"
            ))?
        )
    }
}


// ===============================================================================================
// Python wrapper for Rayleigh process.
// ===============================================================================================
#[pyclass(name = "RayleighProcess", module = "goupil")]
pub struct PyRayleighProcess (RayleighSampler);

#[pymethods]
impl PyRayleighProcess {
    #[new]
    fn new() -> Self {
        Self(RayleighSampler::new(RayleighMode::FormFactor))
    }

    fn cross_section(
        &self,
        energy: Float, // XXX Vectorize these functions.
        material: PyRef<PyMaterialRecord>,
    ) -> Result<Float> {
        let py = material.py();
        let cs = match material.get(py)?.rayleigh_cross_section() {
            None => 0.0,
            Some(table) => table.interpolate(energy),
        };
        Ok(cs)
    }

    fn dcs(
        &self,
        energy: Float,
        cos_theta: Float,
        material: PyRef<PyMaterialRecord>
    ) -> Result<Float> {
        let py = material.py();
        let material = material.get(py)?;
        self.0.dcs(energy, cos_theta, material)
    }

    fn sample(
        &self,
        energy: Float,
        material: PyRef<PyMaterialRecord>
    )
    -> Result<Float> {
        let py = material.py();
        let mut rng = rand::thread_rng();
        let cos_theta = self.0.sample_angle(
            &mut rng,
            energy,
            material.get(py)?
        )?;
        Ok(cos_theta)
    }
}
