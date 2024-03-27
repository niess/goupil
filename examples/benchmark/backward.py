#! /usr/bin/env python3
"""
This example demonstrates the efficiency of backward Monte Carlo methods in
transporting gamma-rays through air. The benchmark experiment involves
calculating the rate of gamma photons collected by a box located 5cm above a
flat ground. It is assumed that radioactive sources are uniformly distributed
throughout the air volume, excluding the collection box.
"""
#! /usr/bin/env python3
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
Note that the WORLD_SIZE parameter could be further increased in a backward
simulation, without loss in Monte Carlo efficiency. However, in order to
reproduce forward simulation results we use the exact same geometry.


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
intentional to achieve a decent Monte Carlo efficiency in the forward case. In
contrast, the efficiency of the backward simulation depends little on the
collector size.

In a backward simulation, gamma-rays start from the collector outer surface and
stop upon reaching a potential source. To maintain consistency with a forward
simulation, backwards transported gamma-rays should be discarded whenever
returning to the collector. This is achieved by setting the box as a transport
boundary, as follows
"""

engine.boundary = collector


"""
================================================================================

  4. Initialising the Monte Carlo states

================================================================================

The last preparatory step is to initialise the Monte Carlo states (representing
gamma-rays). First, we create an empty container of states, as:
"""

N = 1000000
states = goupil.states(N)

"""
In a backward simulation, we specify the final states of photons instead of the
initial ones. This is achieved by using a priori distributions and weighting
them accordingly by the inverse of their PDF.

There are two cases to consider: photo-peaks events, where the final energy
corresponds to the source one, and background events. The set of Monte Carlo
states is randomly split into these two cases according to an a priori fraction
of background events. The DiscreteSpectrum object let us automate these steps,
as
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
source_energies = spectrum.sample(states, engine)

"""
Note that in addition to randomising the final energies, the *sample* function
returns the corresponding source energies.

Arrival positions and directions are randomised over the outer surface of the
collector, using priors. This can be done as
"""

collector.randomise(states, engine)


"""
================================================================================

  5. Running the Monte Carlo simulation

================================================================================

The backward Monte Carlo simulation involves transporting gamma-rays through the
geometry to a potential source. Therefore, we need to specify the anticipated
source energies that establish the stopping criteria for each event.
"""

status = engine.transport(states, source_energies)


"""
================================================================================

  6. Analysing results

================================================================================

The previous transport routines returned an array of status flags indicating the
stop condition for each Monte Carlo state. To estimate the rate of gamma-rays,
we select events that are consistent with a volume source located in the air.
This is done as
"""

sector = geometry.locate(states)
air_index = geometry.sector_index("Atmosphere")
selection = (status == goupil.TransportStatus.ENERGY_CONSTRAINT) & \
            (sector == air_index)
collected = states[selection]

"""
Note that events that exit the simulation geometry or turn back to the collector
are tagged as TransportStatus.EXIT or TransportStatus.BOUNDARY, respectively.

Then, the collected rate is simply obtained by multiplying the Monte Carlo
weights of the selected events with the corresponding activities of sources, per
unit of volume and solid angle, and summing them up. Thus, assuming a volume
activity of 10 Bq/m3 (i.e. 10^5 Bq/cm3 in Goupil's system of units), the source
density is
"""

source_density = 1E-05 / (4.0 * numpy.pi)

"""
where the factor 4 pi accounts for the entire solid angle, assuming an isotropic
emission.
"""

rates = collected["weight"] * source_density / N * 1E-06 # MHz
rate = sum(rates)
sigma_rate = sum(rates**2 - (rate / N)**2)**0.5

"""
Finally, let us print out these results. Note that unlike the forward case, the
rate of collected events and the Monte Carlo efficiency are no longer
proportional.
"""

efficiency = collected.size / N
sigma_efficiency = (efficiency * (1.0 - efficiency) / N)**0.5

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} MHz")
print(f"efficiency = {efficiency:.1E} +- {sigma_efficiency:.1E}")
