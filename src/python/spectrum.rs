use anyhow::Result;
use crate::numerics::{Float, FloatRng};
use crate::transport::{PhotonState, TransportMode};
use pyo3::prelude::*;
use super::macros::value_error;
use super::numpy::{PyArray, PyArrayFlags, PyArrayMethods};
use super::states::States;
use super::transport::PyTransportEngine;
use super::macros::type_error;
use super::rand::PyRandomStream;


// ===============================================================================================
// Source spectrum object.
// ===============================================================================================

#[pyclass(name = "DiscreteSpectrum", module = "goupil")]
pub struct PyDiscreteSpectrum {
    energies: Vec<Float>,
    intensities: Vec<Float>,

    // Backward specific settings.
    #[pyo3(get,set)]
    background: Float,
    #[pyo3(get,set)]
    energy_min: Float,
}

#[pymethods]
impl PyDiscreteSpectrum {
    const DEFAULT_BACKROUND: Float = 0.5;
    const DEFAULT_ENERGY_MIN: Float = 1E-02;

    #[new]
    #[pyo3(signature = (energies, intensities=None, *, background=None, energy_min=None))]
    fn new(
        energies: Vec<Float>,
        intensities: Option<Vec<Float>>,
        background: Option<Float>,
        energy_min: Option<Float>,
    ) -> Result<Self> {
        let intensities = match intensities {
            None => vec![0.0; energies.len()],
            Some(intensities) => if intensities.len() == energies.len() {
                intensities
            } else {
                value_error!(
                    "conflicting energies and intensities sizes ({} != {})",
                    energies.len(),
                    intensities.len(),
                )
            },
        };
        let background = background.unwrap_or(Self::DEFAULT_BACKROUND);
        let energy_min = energy_min.unwrap_or(Self::DEFAULT_ENERGY_MIN);
        let result = Self { energies, intensities, background, energy_min };
        Ok(result)
    }

    #[getter]
    fn get_energies(owner: &Bound<Self>) -> Result<PyObject> {
        let slf = owner.borrow();
        let array = PyArray::from_data(
            owner.py(),
            &slf.energies,
            owner,
            PyArrayFlags::ReadOnly,
            None
        )?;
        Ok(array.into_any().unbind())
    }

    #[getter]
    fn get_intensities(owner: &Bound<Self>) -> Result<PyObject> {
        let slf = owner.borrow();
        let array = PyArray::from_data(
            owner.py(),
            &slf.intensities,
            owner,
            PyArrayFlags::ReadOnly,
            None
        )?;
        Ok(array.into_any().unbind())
    }

    #[pyo3(signature = (states, /, *, engine=None, rng=None, mode=None))]
    fn sample(
        &self,
        states: &Bound<PyAny>,
        engine: Option<&PyTransportEngine>,
        rng: Option<Py<PyRandomStream>>,
        mode: Option<&str>,
    ) -> Result<PyObject> {
        // Unpack arguments.
        let py = states.py();
        let states = States::new(states)?;
        let default_rng: Py<PyRandomStream>;
        let rng = match rng.as_ref() {
            None => match engine.as_ref() {
                None => {
                    default_rng = Py::new(py, PyRandomStream::new(None, None)?)?;
                    &default_rng
                },
                Some(engine) => &engine.random,
            },
            Some(rng) => rng,
        };
        let rng: &mut PyRandomStream = &mut rng.borrow_mut(py);
        let mode = match mode {
            None => match engine.as_ref() {
                None => TransportMode::Forward,
                Some(engine) => {
                    let settings = &engine.settings.borrow(py).inner;
                    settings.mode
                },
            },
            Some(mode) => {
                let mode: TransportMode = mode.try_into()?;
                mode
            },
        };
        if (mode == TransportMode::Backward) && !states.has_weights() {
            type_error!("bad states (expected 'weight' field, found None)")
        }

        // Sample energies.
        let cdf = {
            let mut sum = 0.0;
            let cdf: Vec<Float> = self.intensities
                .iter()
                .map(|v| { sum += v; sum })
                .collect();
            cdf
        };

        let n = states.size();
        let result = match mode {
            TransportMode::Forward => None,
            TransportMode::Backward => {
                let result = PyArray::<Float>::empty(py, &[n])?;
                Some(result)
            },
        };

        for i in 0..n {
            let mut state = states.get(i)?;
            match mode {
                TransportMode::Forward => self.sample_forward(&cdf, rng, &mut state),
                TransportMode::Backward => {
                    let source_energy = self.sample_backward(&cdf, rng, &mut state);
                    result
                        .as_ref()
                        .unwrap()
                        .set(i, source_energy)?;
                },
            }
            states.set(i, &state)?;
        }

        let result = result
            .map(|result| result.into_any().unbind())
            .unwrap_or_else(|| py.None());
        Ok(result)
    }
}

// Private interface.
impl PyDiscreteSpectrum {
    fn sample_backward<R: FloatRng>(
        &self,
        cdf: &[Float],
        rng: &mut R,
        state: &mut PhotonState,
    ) -> Float {
        let line = self.sample_line(cdf, rng);
        if rng.uniform01() < self.background {
            let energy_max = self.energies[line];
            let lne = (energy_max / self.energy_min).ln();
            let energy = self.energy_min * (lne * rng.uniform01()).exp();
            state.energy = energy;
            state.weight *= lne * energy / self.background;
        } else {
            state.energy = self.energies[line];
            state.weight /= 1.0 - self.background;
        }
        self.energies[line]
    }

    fn sample_forward<R: FloatRng>(
        &self,
        cdf: &[Float],
        rng: &mut R,
        state: &mut PhotonState,
    ) {
        let line = self.sample_line(cdf, rng);
        state.energy = self.energies[line];
    }

    fn sample_line<R: FloatRng>(&self, cdf: &[Float], rng: &mut R) -> usize {
        let n = cdf.len();
        let r = cdf[n - 1] * rng.uniform01();
        let line = {
            let mut i = 0;
            loop {
                if r <= cdf[i] { break i }
                else if i == n - 2 { break n - 1 }
                else { i += 1 }
            }
        };
        line
    }
}
