# Using Goupil for a transport problem

This example illustrates the use of Goupil to solve a transport problem using
the Monte Carlo technique. The problem involves calculating the rate of gamma
photons collected by a horizontal surface located 100m above a flat ground. It
is assumed that radioactive sources are uniformly distributed throughout the
ground volume.

In the forward case, the problem's horizontal symmetry is utilized to estimate
the rate of up-going photons. Note that this estimation does not include photons
that would pass the collection surface before being collected and then turn
back. In the backward case, these photons can be easily simulated. However, to
reproduce forward results, the geometry is bounded above 100m by default.

_Note that this example employs a [stratified geometry][STRATIFIED_GEOMETRY]
with a flat ground interface._


## Usage

The Python scripts [forward.py](forward.py) and [backward.py](backward.py)
execute a Monte Carlo simulation in either the forward or backward direction. In
the backward case, if the collection surface is removed from the geometry, then
photons that turn back above 100m are simulated.


[STRATIFIED_GEOMETRY]: https://goupil.readthedocs.io/en/latest/py/stratified_geometry.html
