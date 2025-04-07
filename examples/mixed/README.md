# Mixed Monte Carlo simulation

This example illustrates the mixed simulation scheme described in
[NVT24][NVT24]. The test [geometry](geometry.toml) consists of a cylindrical
scintillation detector immersed in water. The response of the detector to
[radio-isotopes](sources.toml) that are uniformly scattered throughout the water
volume is computed in two steps. Firstly, energy deposits are simulated over the
scintillator volume in response to an incoming gamma-ray flux using Geant4
(through [Calzone][CALZONE]). Then, the incoming flux (and thus energy deposits)
are weighted according to radioactive sources using Goupil in backward Monte
Carlo mode. Finally, the detector energy response is modelled by a Gausian
distribution, applied to the Monte Carlo energy deposits.

## Usage

In order to run this example, it is first necessary to install
[Calzone][CALZONE] (available from [PyPI][PYPI]), a Python Geant4 wrapper that
is interoperable with Goupil. The Python script [run-mixed.py](run-mixed.py) (or
[run-forward.py](run-forward.py)) executes a mixed (end-to-end forward) Monte
Carlo simulation. By default, the energy deposits are recorded under the
`data/mixed` (`data/forward`) folder. The
[compute-response.py](compute-response.py) script is then used to calculate the
detector's response to the input energy deposits. The behaviour of these scripts
is configurable through Command Line Interfaces (CLIs). For instance, the
following commands will run `100,000` Monte Carlo events in mixed mode, compute
the corresponding energy response, and display the result (using matplotlib):

```bash
./run-mixed.py -n 100000 && ./compute-response.py -p
```

[CALZONE]: https://github.com/niess/calzone
[NVT24]: https://doi.org/10.48550/arXiv.2412.02414
[PYPI]: https://pypi.org/project/calzone/
