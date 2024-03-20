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

# Create a stratified geometry covering a 20x20 km^2 area. Rocks and air are
# separated by a plane topography surface at z = 0. In addition, we bound the
# geometry above with a plane collection surface.

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
    collection_surface,
    goupil.GeometrySector(air, air_density, "Atmosphere"),
    topography_surface,
    goupil.GeometrySector(rock, rock_density, "Ground")
)

# Create a transport engine (configured in forward mode, by default).
engine = goupil.TransportEngine(geometry)

# Create an array of Monte Carlo states.
N = 1000000
states = goupil.states(N)

# Randomise initial energies using the main emission lines of Pb-214 decay.
source_spectrum = numpy.array((
    (0.242, 0.12),
    (0.295, 0.30),
    (0.352, 0.58),
))
states["energy"] = numpy.random.choice(
    source_spectrum[:,0],
    N,
    replace = True,
    p = source_spectrum[:,1]
)

# Randomise the source depth.
MAX_DEPTH = 1.0E+02 # cm
states["position"][:,2] = -MAX_DEPTH * engine.random.uniform01(N)

# Randomise the emission direction (uniformly over the entire solid angle).
cos_theta = 2.0 * engine.random.uniform01(N) - 1.0
sin_theta = numpy.sqrt(1.0 - cos_theta**2)
phi = 2.0 * numpy.pi * engine.random.uniform01(N)

states["direction"][:,0] = numpy.cos(phi) * sin_theta
states["direction"][:,1] = numpy.sin(phi) * sin_theta
states["direction"][:,2] = cos_theta

# Run the Monte Carlo transport.
status = engine.transport(states)

# Select upgoing events that exit through the collection surface with an energy
# greater than 10 keV.

ENERGY_MIN = 1E-02 # MeV
selection = (status == goupil.TransportStatus.EXIT) & \
            (states["position"][:,2] >= COLLECTION_HEIGHT) & \
            (states["direction"][:,2] > 0.0) & \
            (states["energy"] >= ENERGY_MIN)
collected = states[selection]

# Compute the rate of events crossing the collection surface, assuming a volume
# activity of 1 Bq / cm^2 sr is assumed.

volume_activity = 1.0 # Bq / (cm^3 sr)
source_volume = (2.0 * HALF_WIDTH)**2 * MAX_DEPTH
solid_angle = 4.0 * numpy.pi
total_activity = volume_activity * source_volume * solid_angle

collection_area = (2.0 * HALF_WIDTH)**2

p = collected.size / N
rate = p * total_activity / collection_area

# Print results.
#
# Note that in the forward case, the rates of collected events is directly
# proportional to the Monte Carlo efficiency.

efficiency = 100.0 * p
sigma = numpy.sqrt(p * (1.0 - p) / N)
sigma_efficiency = 100.0 * sigma
sigma_rate = total_activity / collection_area * sigma

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} Bq / cm^2")
print(f"efficiency = {efficiency:.2f} +- {sigma_efficiency:.2f} %")