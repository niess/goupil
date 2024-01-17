use anyhow::Result;
use crate::numerics::Float;
use crate::transport::{
    density::DensityModel,
    geometry::{ExternalGeometry, ExternalTracer, GeometryDefinition, GeometryTracer,
               SimpleGeometry},
    PhotonState,
};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use super::ctrlc_catched;
use super::macros::value_error;
use super::materials::PyMaterialDefinition;
use super::numpy::{ArrayOrFloat, PyArray};
use super::transport::CState;


// ===============================================================================================
// Python wrapper for a simple geometry object.
// ===============================================================================================

#[pyclass(name = "SimpleGeometry", module = "goupil")]
pub struct PySimpleGeometry (pub SimpleGeometry);

#[pymethods]
impl PySimpleGeometry {
    #[new]
    fn new(
        material: PyRef<PyMaterialDefinition>,
        density: DensityModel,
    ) -> Result<Self> {
        let geometry = SimpleGeometry::new(&material.0, density);
        Ok(Self(geometry))
    }

    #[getter]
    fn get_density(&self, py: Python) -> PyObject {
        self.0.sectors[0].density.into_py(py)
    }

    #[setter]
    fn set_density(&mut self, value: DensityModel) -> Result<()> {
        self.0.sectors[0].density = value;
        Ok(())
    }

    #[getter]
    fn get_material(&self) -> PyMaterialDefinition {
        PyMaterialDefinition(self.0.materials[0].clone())
    }
}


// ===============================================================================================
// Python wrapper for an external geometry object.
// ===============================================================================================

#[pyclass(name = "ExternalGeometry", module = "goupil")]
pub struct PyExternalGeometry (pub ExternalGeometry);

#[pymethods]
impl PyExternalGeometry {
    #[new]
    pub fn new(path: &str) -> Result<Self> {
        let geometry = unsafe { ExternalGeometry::new(path)? };
        Ok(Self(geometry))
    }

    #[getter]
    fn get_materials<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let mut materials = Vec::<PyObject>::with_capacity(self.0.materials.len());
        for material in self.0.materials.iter() {
            let material = PyMaterialDefinition(material.clone());
            materials.push(material.into_py(py));
        }
        PyTuple::new(py, materials)
    }

    #[getter]
    fn get_sectors<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let sectors: Vec<_> = self.0
            .sectors
            .iter()
            .map(|sector| (
                sector.material,
                sector.density,
                sector.description
                    .as_ref()
                    .map(|description| description.to_string()),
            ))
            .collect();
        PyTuple::new(py, sectors)
    }

    fn locate(&self, states: &PyArray<CState>) -> Result<PyObject> {
        let py = states.py();
        let sectors = PyArray::<usize>::empty(py, &states.shape())?;
        let mut tracer = ExternalTracer::new(&self.0)?;
        let m = self.0.sectors().len();
        let n = states.size();
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            tracer.reset(state.position, state.direction)?;
            let sector = tracer.sector().unwrap_or(m);
            sectors.set(i, sector)?;

            if i % 1000 == 0 { // Check for a Ctrl+C interrupt, catched by Python.
                ctrlc_catched()?;
            }
        }
        let sectors: &PyAny = sectors;
        Ok(sectors.into_py(py))
    }

    fn trace(
        &self,
        states: &PyArray<CState>,
        lengths: Option<ArrayOrFloat>,
        density: Option<bool>,
    ) -> Result<PyObject> {
        let n = states.size();
        if let Some(lengths) = lengths.as_ref() {
            if let ArrayOrFloat::Array(lengths) = &lengths {
                if lengths.size() != states.size() {
                    value_error!(
                        "bad lengths (expected a float or a size {} array, found a size {} array)",
                        states.size(),
                        lengths.size(),
                    )
                }
            }
        }

        let mut shape = states.shape();
        let m = self.0.sectors().len();
        shape.push(m);
        let py = states.py();
        let result = PyArray::<Float>::empty(py, &shape)?;

        let density = density.unwrap_or(false);
        let mut tracer = ExternalTracer::new(&self.0)?;
        let mut k: usize = 0;
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            let mut grammages: Vec<Float> = vec![0.0; m];
            tracer.reset(state.position, state.direction)?;
            let mut length = match lengths.as_ref() {
                None => Float::INFINITY,
                Some(lengths) => match &lengths {
                    ArrayOrFloat::Array(lengths) => lengths.get(i)?,
                    ArrayOrFloat::Float(lengths) => *lengths,
                },
            };
            loop {
                match tracer.sector() {
                    None => break,
                    Some(sector) => {
                        let step_length = tracer.trace(length)?;
                        if density {
                            let model = &self.0.sectors[sector].density;
                            let position = tracer.position();
                            grammages[sector] += model.column_depth(
                                position, state.direction, step_length
                            );
                        } else {
                            grammages[sector] += step_length;
                        }
                        if lengths.is_some() {
                            length -= step_length;
                            if length <= 0.0 { break }
                        }
                        tracer.update(step_length, state.direction)?;
                    },
                }
                k += 1;
                if k == 1000 { // Check for a Ctrl+C interrupt, catched by Python.
                    ctrlc_catched()?;
                    k = 0;
                }
            }
            let j0 = i * m;
            for j in 0..m {
                result.set(j0 + j, grammages[j])?;
            }
        }
        let result: &PyAny = result;
        Ok(result.into_py(py))
    }

    fn update_material(
        &mut self,
        index: usize,
        material: PyRef<PyMaterialDefinition>
    ) -> Result<()> {
        self.0.update_material(index, &material.0)
    }

    fn update_sector(
        &mut self,
        index: usize,
        material: Option<usize>,
        density: Option<DensityModel>,
    ) -> Result<()> {
        self.0.update_sector(index, material, density.as_ref())
    }
}


// ===============================================================================================
// Unresolved geometry definition.
// ===============================================================================================

#[derive(Clone, FromPyObject)]
pub enum PyGeometryDefinition {
    External(Py<PyExternalGeometry>),
    Simple(Py<PySimpleGeometry>),
}

impl IntoPy<PyObject> for PyGeometryDefinition {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::External(external) => external.into_py(py),
            Self::Simple(simple) => simple.into_py(py),
        }
    }
}
