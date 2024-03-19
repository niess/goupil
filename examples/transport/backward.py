#! /usr/bin/env python3
import goupil
import numpy


# Define the air material.
air = goupil.MaterialDefinition(
    "Air",
    mass_composition = (
        (0.76, "N"),
        (0.23, "O"),
        (0.01, "Ar"),
    )
)
air_density = goupil.DensityGradient(
    1.225E-03, # g/cm^3
    1.04E+06   # cm
)

# Define the rocks material.
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

# Create a stratified geometry covering a 2x2 km^2 area. Rocks and air are
# separated by a plane topography surface at z = 0. In addition, we bound the
# geometry above with a plane collection surface.

HALF_WIDTH, COLLECTION_HEIGHT = 1E+05, 1E+03 # cm
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
    collection_surface,
    goupil.GeometrySector(air, air_density, "Atmosphere"),
    topography_surface,
    goupil.GeometrySector(rock, rock_density, "Ground")
)

# Create a transport engine, configured in backward mode.
engine = goupil.TransportEngine(geometry)
engine.mode = "Backward"

# Create an array of Monte Carlo states.
N = 1000000
states = goupil.states(N)

# Randomise initial energies using the main emission lines of Pb-214 decay.
source_spectrum = numpy.array((
    (0.242, 0.12),
    (0.295, 0.30),
    (0.352, 0.58),
))
source_energies = numpy.random.choice(
    source_spectrum[:,0],
    N,
    replace = True,
    p = source_spectrum[:,1]
)

ALPHA = 0.5 # Fraction of photo-peaks events.
photopeaks = engine.random.uniform01(N) < ALPHA
states["energy"][photopeaks] = source_energies[photopeaks]
states["weight"][photopeaks] = 1.0 / ALPHA

background = ~photopeaks
ENERGY_MIN = 1E-02 # MeV
lne = numpy.log(source_energies[background] / ENERGY_MIN)
energies = ENERGY_MIN * numpy.exp(lne * engine.random.uniform01(lne.size))
states["energy"][background] = energies
states["weight"][background] = lne * energies / (1.0 - ALPHA)

# Randomise the arrival position.
epsilon = 1E-04 # cm, for numerical safety.
states["position"][:,0] = HALF_WIDTH * (2.0 * engine.random.uniform01(N) - 1.0)
states["position"][:,1] = HALF_WIDTH * (2.0 * engine.random.uniform01(N) - 1.0)
states["position"][:,2] = COLLECTION_HEIGHT - epsilon
states["weight"] *= (2.0 * HALF_WIDTH)**2

# Randomise the arrival direction (upgoing).
#
# A cosine distribution is used in order to account for the crossing-factor
# (cos(theta)) of the collection surface (see below).

u = engine.random.uniform01(N)
cos_theta = numpy.sqrt(u)
sin_theta = numpy.sqrt(1.0 - u)
phi = 2.0 * numpy.pi * engine.random.uniform01(N)

states["direction"][:,0] = numpy.cos(phi) * sin_theta
states["direction"][:,1] = numpy.sin(phi) * sin_theta
states["direction"][:,2] = cos_theta

states["weight"] *= numpy.pi # Note that the angular weight includes a surface
                             # crossing factor (cos(theta)) which results from
                             # the flux definition. This factor simplifies out
                             # with the generation PDF (cos(theta) / pi).

# Run the Monte Carlo transport.
status = engine.transport(
    states,
    source_energies
)

# Select events that are consistent with a source located in the ground, within
# 1m depth.

MAX_DEPTH = 1.0E+02 # cm
sector = geometry.locate(states)
sector_names = [sector.description for sector in geometry.sectors]
ground_index = sector_names.index("Ground")
selection = (status == goupil.TransportStatus.ENERGY_CONSTRAINT) & \
            (sector == ground_index) & \
            (states["position"][:,2] >= -MAX_DEPTH)
collected = states[selection]

# Print statistics.
source_volume = (2.0 * HALF_WIDTH)**2 * MAX_DEPTH
source_density = 1.0 / (4.0 * numpy.pi * source_volume) # A normalised source
                                                        # intensity is assumed.

rates = collected["weight"] * source_density / N
mu = 100.0 * sum(rates)
sigma = 100.0 * sum(rates**2 - (sum(rates) / N)**2)**0.5
print(f"rate = {mu:.2f} +- {sigma:.2f} %")

m = collected.size
efficiency = 100.0 * m / N
sigma = 100.0 * numpy.sqrt(m * (1.0 - m / N)) / N
print(f"efficiency = {efficiency:.2f} +- {sigma:.2f} %")
