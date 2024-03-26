#! /usr/bin/env python3
"""
This example illustrates the use of Goupil to solve a transport problem using
the backward Monte Carlo technique. The problem involves calculating the rate of
gamma photons collected by a horizontal surface located 100m above a flat
ground. It is assumed that radioactive sources are uniformly distributed
throughout the ground volume.
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
topography surface at z = 0. Additionally, the geometry is bounded above by a
plane representing the collection surface for gamma-rays.

Let us point out that the collection surface is not necessary for a backward
simulation. This boundary prevents the simulation of photons that would pass the
collection surface, turn back, and finally be collected as up-going. To maintain
consistency with the forward Monte Carlo simulation, these events are not
included. However, if you wish to simulate them, simply remove the collection
surface boundary.
"""

HALF_WIDTH, COLLECTION_HEIGHT = 1E+06, 1E+03 # cm
topography_surface = goupil.TopographyMap(
    (-HALF_WIDTH, HALF_WIDTH),
    (-HALF_WIDTH, HALF_WIDTH),
    z = 0 # cm
)
collection_surface = goupil.TopographyMap(
    (-HALF_WIDTH, HALF_WIDTH),
    (-HALF_WIDTH, HALF_WIDTH),
    z = COLLECTION_HEIGHT
)
geometry = goupil.StratifiedGeometry(
    collection_surface, # <== Comment this line in order to include photons that
                        #     would go past the collection height, turn back,
                        #     and finally be collected as up-going.
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
engine.mode = "Backward"

"""
Note that the transport engine is configured for forward Monte Carlo by default.
So we needed to adjust its mode of operation. Note also that the engine provides
a random stream which we will use below to initialise Monte Carlo states.
"""


"""
================================================================================

  3. Initialising the Monte Carlo states

================================================================================

The second preparatory step is to initialise the Monte Carlo states
(representing gamma-rays). First, we create an empty container of states, as:
"""

N = 1000000
states = goupil.states(N)

"""
For this problem we consider a single radioactive isotope, Pb-214, with its
three main gamma emission lines. The corresponding spectrum is
"""

source_spectrum = numpy.array((
    (0.242, 0.12),
    (0.295, 0.30),
    (0.352, 0.58),
))

"""
Let us randomise the source energies of the gamma-rays according to the latter
spectrum. This is done using the `random_choice` function from the numpy
package.
"""

source_energies = numpy.random.choice(
    source_spectrum[:,0],
    N,
    replace = True,
    p = source_spectrum[:,1]
)

"""
In a backward simulation, we specify the final states of photons instead of the
initial ones. This is achieved by using a priori distributions and weighting
them accordingly by the inverse of their PDF.

There are two cases to consider: photo-peaks events, where the final energy
corresponds to the source one, and background events. The set of Monte Carlo
states is randomly split into these two cases with a probability alpha. Setting
the energy of photo-peaks is straightforward, as follows
"""

ALPHA = 0.5 # Fraction of photo-peaks events.
photopeaks = engine.random.uniform01(N) < ALPHA
states["energy"][photopeaks] = source_energies[photopeaks]
states["weight"][photopeaks] = 1.0 / ALPHA # i.e., 1 / PDF

"""
For background events, we use a log-uniform prior to set the final energy.
"""

background = ~photopeaks
ENERGY_MIN = 1E-02 # MeV
lne = numpy.log(source_energies[background] / ENERGY_MIN)
energies = ENERGY_MIN * numpy.exp(lne * engine.random.uniform01(lne.size))
states["energy"][background] = energies
states["weight"][background] = lne * energies / (1.0 - ALPHA) # i.e., 1 / PDF

"""
Arrival positions are constant and thus easily set. To avoid numerical problems,
a small offset is applied with respect to the nominal position, which would be
on the collection surface.
"""

epsilon = 1E-04 # cm, for numerical safety.
states["position"][:,0] = 0.0
states["position"][:,1] = 0.0
states["position"][:,2] = COLLECTION_HEIGHT - epsilon

"""
The directions of arrival are randomised over the upper half of the solid angle
because we only collect up-going photons.

A cosine distribution is used to account for the crossing-factor (cos(theta)) of
the collection surface. That is, the angular weight below includes the cosine
factor related to the intersection of the collection surface, not just the PDF
prior.
"""

u = engine.random.uniform01(N)
cos_theta = numpy.sqrt(u)
sin_theta = numpy.sqrt(1.0 - u)
phi = 2.0 * numpy.pi * engine.random.uniform01(N)

states["direction"][:,0] = numpy.cos(phi) * sin_theta
states["direction"][:,1] = numpy.sin(phi) * sin_theta
states["direction"][:,2] = cos_theta

states["weight"] *= numpy.pi # Note that the angular weight includes a surface
                             # crossing factor (cos(theta)). This factor
                             # simplifies out with the generation PDF
                             # (cos(theta) / pi).


"""
================================================================================

  4. Running the Monte Carlo simulation

================================================================================

The backward Monte Carlo simulation involves transporting gamma-rays through the
geometry to a potential source. Therefore, we need to specify the anticipated
source energies that establish the stopping criteria for each event.
"""

status = engine.transport(
    states,
    source_energies
)


"""
================================================================================

  5. Analysing results

================================================================================

The previous transport routines returned an array of status flags indicating the
stop condition for each Monte Carlo state. To estimate the rate of gamma-rays,
we select events that are consistent with a source located below the ground,
within a depth of 1m.
"""

MAX_DEPTH = 1.0E+02 # cm
sector = geometry.locate(states)
ground_index = geometry.sector_index("Ground")
selection = (status == goupil.TransportStatus.ENERGY_CONSTRAINT) & \
            (sector == ground_index) & \
            (states["position"][:,2] >= -MAX_DEPTH)
collected = states[selection]

"""
Then, the collected rate is simply obtained by multiplying the Monte Carlo
weights of the selected events with the corresponding volume activities of
sources, and summing them up. Thus
"""

volume_activity = 1.0 # Bq / (cm^3 sr)

rates = collected["weight"] * volume_activity / N
rate = sum(rates)
sigma_rate = ((sum(rates**2) - rate**2 / N))**0.5

"""
Finally, let us print out these results. Note that unlike the forward case, the
rate of collected events and the Monte Carlo efficiency are no longer
proportional.
"""

p = collected.size / N
efficiency = 100.0 * p
sigma_efficiency = 100.0 * numpy.sqrt(p * (1.0 - p) / N)

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} Hz / cm^2")
print(f"efficiency = {efficiency:.2f} +- {sigma_efficiency:.2f} %")
