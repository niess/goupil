#! /usr/bin/env python3
"""
This example illustrates the use of Goupil to solve a transport problem using
the conventional (forward) Monte Carlo technique. The problem involves
calculating the rate of gamma photons collected by a horizontal surface located
100m above a flat ground. It is assumed that radioactive sources are uniformly
distributed throughout the ground volume.
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
Secondly, we define the composition of the rocks that make up the lower volume
of the Monte Carlo geometry. In this case, a uniform density is assumed.
"""

rock = goupil.MaterialDefinition(
    "Rock",
    mass_composition = (
        (0.57, "SiO2"),
        (0.19, "Al2O3"),
        (0.12, "FeO"),
        (0.12, "CaO"),
    )
)
rock_density = 2.9 # g/cm^3

"""
Finally, a vertically stratified Monte Carlo geometry is created with a 20x20
km^2 horizontal base. The rocks and air volumes are separated by a plane
topography surface at z = 0.
"""

WORLD_WIDTH = 2E+06 # cm
topography_surface = goupil.TopographyMap(
    (-0.5 * WORLD_WIDTH, 0.5 * WORLD_WIDTH),
    (-0.5 * WORLD_WIDTH, 0.5 * WORLD_WIDTH),
    z = 0 # cm
)
geometry = goupil.StratifiedGeometry(
    goupil.GeometrySector(air, air_density, "Atmosphere"),
    topography_surface,
    goupil.GeometrySector(rock, rock_density, "Ground")
)


"""
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
"""


"""
================================================================================

  3. Initialising the Monte Carlo states

================================================================================

The second preparatory step is to initialise the Monte Carlo states
(representing gamma-rays) according to the distribution of radioactive sources.
First, we create an empty container of states, as:
"""

N = 1000000
states = goupil.states(N)

"""
For this problem we consider a single radioactive isotope, Pb-214, with its
three main gamma emission lines. The corresponding spectrum is
"""

spectrum = goupil.DiscreteSpectrum(
    energies = (0.242, 0.295, 0.352), # MeV
    intensities = (7.3, 18.4, 35.6)
)

"""
Let us randomise the initial energies of the gamma-rays according to the latter
spectrum. This is done as follows
"""

spectrum.sample(states, engine)

"""
In principle, the source positions should be randomized throughout the entire
ground volume and counted at (0, 0, h), where h represents the collection
height. However, due to the problem's horizontal symmetry, it is equivalent to
emit photons from (0, 0, z0), where z0 < 0, and count them over the entire
collection surface. This property is exploited to enhance the efficiency of the
Monte Carlo. Therefore, sources are only generated along a vertical column
centered on the origin, horizontally.
"""

MAX_DEPTH = 1.0E+02 # cm
states["position"][:,2] = -MAX_DEPTH * engine.random.uniform01(N)

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

  4. Running the Monte Carlo simulation

================================================================================

The Monte Carlo simulation is carried out by transporting the emitted gamma-rays
through the geometry. But, first we need to define a boundary plane for
collecting gamma rays. This is done with a 1mm thin box, as follows:
"""

COLLECTION_HEIGHT, COLLECTION_THICKNESS = 1E+03, 0.1  # cm
engine.boundary = goupil.BoxShape(
    size = (WORLD_WIDTH, WORLD_WIDTH, COLLECTION_THICKNESS), # cm
    center = (0.0, 0.0, COLLECTION_HEIGHT + 0.5 * COLLECTION_THICKNESS) # cm
)

"""
Then, the Monte Carlo runs as
"""

status = engine.transport(states)


"""
================================================================================

  5. Analysing results

================================================================================

The previous transport routines returned an array of status flags indicating the
stop condition for each Monte Carlo state. To estimate the rate of gamma-rays at
point (0,0,h), we select up-going photons that entered the collection
surface with an energy higher than 10 keV.
"""

ENERGY_MIN = 1E-02 # MeV
selection = (status == goupil.TransportStatus.BOUNDARY) & \
            (states["direction"][:,2] > 0.0) & \
            (states["energy"] >= ENERGY_MIN)
collected = states[selection]

"""
The number of collected photons provides a Monte Carlo estimate of the rate of
gamma-rays. To compute this rate, we need to determine the total activity of
sources represented by our Monte Carlo simulation. Assuming a unit volume
activity, the total activity is obtained by integrating over the source volume
and solid angle.
"""

volume_activity = 1.0 # Bq / (cm^3 sr)
source_volume = WORLD_WIDTH**2 * MAX_DEPTH
solid_angle = 4.0 * numpy.pi
total_activity = volume_activity * source_volume * solid_angle

"""
The activity registered on the collection surface is a fraction k/N of the total
source activity, where k is the number of collected photons. The rate of gamma
rays is then calculated by dividing the number of collected photons by the
collection area. Therefore, we define a rate factor
"""

collection_area = WORLD_WIDTH**2
K = total_activity / collection_area

"""
from which we obtain the rate of collected photons
"""

p = collected.size / N
rate = p * K

"""
Finally, let us print out these results. Note that in the forward case, the
rates of collected events and the Monte Carlo efficiency are proportional.
"""

efficiency = 100.0 * p
sigma = numpy.sqrt(p * (1.0 - p) / N)
sigma_efficiency = 100.0 * sigma
sigma_rate = K * sigma

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} Hz / cm^2")
print(f"efficiency = {efficiency:.2f} +- {sigma_efficiency:.2f} %")
