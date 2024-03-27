# Using Goupil with a Geant4 geometry

This example demonstrates how Geant4 can be used as an external geometry engine
for Goupil. The geometry is exported as a shared library and then loaded by
Goupil as an [ExternalGeometry][EXTERNAL_GEOMETRY] Python object. In addition,
the geometry library includes custom functions for initializing Monte Carlo
states. These functions are accessed from Python by using the [ctypes][CTYPES]
module.

_Note that this example uses the [Geant4 adapter][G4_GOUPIL] distributed with
Goupil. An alternative implementation of this problem (not requiring Geant4) can
be found under [examples/benchmark](../benchmark)._


## Usage

To run this example, you must first compile the geometry library
(_libgeometry.so_), for instance, by using the provided [Makefile](Makefile)
(assuming a Unix system). The Python script [forward.py](forward.py) (or
[backward.py](backward.py)) executes a Monte Carlo simulation in the
forward (or backward) direction, using the compiled Geant4 geometry.

_Note that compiling and running the geometry library requires both Geant4
and Goupil._


[CTYPES]: https://docs.python.org/3/library/ctypes.html
[EXTERNAL_GEOMETRY]: https://goupil.readthedocs.io/en/latest/py/external_geometry.html
[G4_GOUPIL]: https://goupil.readthedocs.io/en/latest/geant4.html
