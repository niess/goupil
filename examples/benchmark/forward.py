#! /usr/bin/env python3
"""
This example demonstrates the inefficiency of conventional (analogue) Monte
Carlo methods in transporting gamma-rays through air. The benchmark experiment
involves calculating the rate of gamma photons collected by a box located 5cm
above a flat ground. It is assumed that radioactive sources are uniformly
distributed throughout the air volume, excluding the collection box.
"""
import goupil
import numpy


"""
================================================================================

  1. Defining the Monte Carlo geometry

================================================================================

To begin, we define the geometry of the Monte Carlo simulation. Firstly, we
define the material that makes up the upper air volume.
"""

air = goupil.MaterialDefinition(
    "Air",
    mass_composition = (
        (0.76, "N"),
        (0.23, "O"),
        (0.01, "Ar"),
    )
)

"""
The air density is modelled by a vertical gradient, decreasing exponentially
with altitude. This is a suitable approximation for the lower layer of the
atmosphere, i.e. the troposphere.
"""

air_density = goupil.DensityGradient(
    1.225E-03, # g/cm^3
    1.04E+06   # cm
)

"""
Secondly, we define the composition of the soil that makes up the lower volume
of the Monte Carlo geometry. In this case, a Limestone soil of uniform density
is assumed.
"""

rock = "CaCO3"
rock_density = 2.8 # g/cm^3

"""
Finally, a vertically stratified Monte Carlo geometry is created with a 2x2
km^2 horizontal base. The rocks and air volumes are separated by a plane
topography surface at z = 0.
"""

WORLD_SIZE = 2E+05 # cm
interface = goupil.TopographyMap(
    (-0.5 * WORLD_SIZE, 0.5 * WORLD_SIZE),
    (-0.5 * WORLD_SIZE, 0.5 * WORLD_SIZE),
    z = 0 # cm
)
geometry = goupil.StratifiedGeometry(
    goupil.GeometrySector(air, air_density, "Atmosphere"),
    interface,
    goupil.GeometrySector(rock, rock_density, "Ground")
)

"""
It should be noted that the WORLD_SIZE parameter has a direct impact on the
Monte Carlo efficiency of the forward simulation. To achieve a 1% accuracy on
the collection rate, an extent of +- 1km around the detector box is required.
Decreasing the world size leads to a significant loss (due to some remote gamma
rays being ignored), while increasing it further results in an efficiency drop
without a significant increase in the collection rate.


================================================================================

  2. Preparing the transport engine

================================================================================

The Monte Carlo simulation is performed by a transport engine. Let us create
this engine and initialise it according to the previous Monte Carlo geometry.
This is done as follows.
"""

engine = goupil.TransportEngine(geometry)

"""
Note that the transport engine is configured for forward Monte Carlo by default.
So we do not need to adjust it further. Note also that the engine provides a
random stream which we will use below to initialise Monte Carlo states.


================================================================================

  3. Setting the collection box

================================================================================

The collector is defined as a 20x20x10 m3 box. This is done as
"""

DETECTOR_WIDTH, DETECTOR_HEIGHT, DETECTOR_OFFSET = 2E+03, 1E+03, 5.0 # cm
collector = goupil.BoxShape(
    size = (DETECTOR_WIDTH, DETECTOR_WIDTH, DETECTOR_HEIGHT),
    center = (0.0, 0.0, 0.5 * DETECTOR_HEIGHT + DETECTOR_OFFSET)
)

"""
Note that the collector is larger than a typical gamma-ray spectrometer. This is
intentional to achieve a decent Monte Carlo efficiency in the forward case.

Gamma-rays are stopped upon entering the collector by setting the box as a
transport boundary, as follows
"""

engine.boundary = collector


"""
================================================================================

  4. Initialising the Monte Carlo states

================================================================================

The last preparatory step is to initialise the Monte Carlo states
(representing gamma-rays) according to the distribution of radioactive sources.
First, we create an empty container of states, as:
"""

N = 1000000
states = goupil.states(N)

"""
For this problem we consider two radioactive isotopes, Pb-214 and Bi-214 (i.e.
Radon-222 progenies), with their main gamma emission lines. The corresponding
spectrum is
"""

spectrum = numpy.array((
    # Pb-214 major emission lines.
    (0.242,  7.3),
    (0.295, 18.4),
    (0.352, 35.6),
    # Bi-214 major emission lines.
    (0.609, 45.5),
    (0.768,  4.9),
    (0.934,  3.1),
    (1.120, 14.9),
    (1.238,  5.8),
    (1.378,  4.0),
    (1.764, 15.3),
    (2.204,  4.9),
))
spectrum = goupil.DiscreteSpectrum(
    energies = spectrum[:,0], # MeV
    intensities = spectrum[:,1],
)

"""
Let us randomise the initial energies of the gamma-rays according to the latter
spectrum. This is done as follows
"""

spectrum.sample(states, engine=engine)

"""
The positions of the sources are uniformly randomised throughout the air volume,
including the collector volume at this stage.
"""

states["position"][:,0] = WORLD_SIZE * (engine.random.uniform01(N) - 0.5)
states["position"][:,1] = WORLD_SIZE * (engine.random.uniform01(N) - 0.5)
states["position"][:,2] = 0.5 * WORLD_SIZE * engine.random.uniform01(N)

"""
Next, sources inside the collector are discarded. This only applies to a small
number of events, as the world volume is significantly larger than the detector
volume.
"""

states = states[~collector.inside(states)]
N = states.size

"""
Gamma-rays are assumed to be emitted isotropically. Thus, their initial
directions are generated uniformly over the entire solid angle.
"""

cos_theta = 2.0 * engine.random.uniform01(N) - 1.0
sin_theta = numpy.sqrt(1.0 - cos_theta**2)
phi = 2.0 * numpy.pi * engine.random.uniform01(N)

states["direction"][:,0] = numpy.cos(phi) * sin_theta
states["direction"][:,1] = numpy.sin(phi) * sin_theta
states["direction"][:,2] = cos_theta


"""
================================================================================

  5. Running the Monte Carlo simulation

================================================================================

The Monte Carlo simulation is carried out by transporting the emitted gamma-rays
through the geometry. This is simply done as follows
"""

status = engine.transport(states)


"""
================================================================================

  6. Analysing results

================================================================================

The previous transport routines returned an array of status flags indicating the
stop condition for each Monte Carlo state. To estimate the rate of collected
gamma-rays we select events that reached the box boundary, as
"""
collected = states[status == goupil.TransportStatus.BOUNDARY]

"""
Note that events that exit the simulation geometry are tagged as
TransportStatus.EXIT instead of TransportStatus.BOUNDARY. The latter flag only
applies to the explicit engine.boundary condition.

The number of collected photons provides a Monte Carlo estimate of the rate of
gamma-rays. To compute this rate, we need to determine the total activity of
sources represented by our Monte Carlo simulation. Assuming a volume activity of
10 Bq/m3 (i.e. 10^5 Bq/cm3 in Goupil's system of units), the total activity is
"""

source_density = 1E-05 # Bq/cm^3
source_volume = 0.5 * WORLD_SIZE**3 - DETECTOR_WIDTH**2 * DETECTOR_HEIGHT
total_activity = source_density * source_volume

"""
The activity registered on the collection surface is a fraction k/N of the total
source activity, where k is the number of collected photons. Thus
"""

efficiency = collected.size / N
sigma_efficiency = (efficiency * (1.0 - efficiency) / N)**0.5

rate = efficiency * total_activity * 1E-06 # MHz
sigma_rate = sigma_efficiency * total_activity * 1E-06 # MHz

"""
Finally, let us print out these results. Note that in the forward case, the
rates of collected events and the Monte Carlo efficiency are proportional.
"""

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} MHz")
print(f"efficiency = {efficiency:.1E} +- {sigma_efficiency:.1E}")
