# Benchmarking Monte Carlo efficiency

This example demonstrates the use of Goupil to enhance the Monte Carlo
efficiency for transporting gamma rays through the air.

The geometry being considered is a flat soil-air interface with a box-like
detector positioned 5cm above the ground. Radioactive sources are uniformly
distributed throughout the air volume surrounding the box. The benchmark
experiment involves counting gamma-rays collected over the surface of the
detector box.

In the forward case, the collection rate is proportional to the Monte Carlo
efficiency, which is very low. However, in the backward case, the collection
rate can be accurately estimated while maintaining a high Monte Carlo
efficiency.

_This example uses a [stratified geometry][STRATIFIED_GEOMETRY]. An alternative
implementation, which utilizes Geant4's ray-tracer, can be found in the
[examples/geant4](../geant4) directory._


## Usage

The Python scripts [forward.py](forward.py) and [backward.py](backward.py)
execute a Monte Carlo simulation in either the forward or backward direction.


[STRATIFIED_GEOMETRY]: https://goupil.readthedocs.io/en/latest/py/stratified_geometry.html
