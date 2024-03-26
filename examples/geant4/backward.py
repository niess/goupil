#! /usr/bin/env python3
import ctypes
import goupil
import numpy


# Load the external geometry (using an explicit relative path).
geometry = goupil.ExternalGeometry("./libgeometry.so")

# Configure the user defined external function(s) (using ctypes).
geometry.lib.initialise_states_backward.argtypes = [
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.c_void_p,
    ctypes.c_void_p,
]
geometry.lib.initialise_states_backward.restype = None

# Set a vertical density gradient for the atmosphere.
geometry.update_sector(
    sector = "Atmosphere",
    density = goupil.DensityGradient(
        1.225E-03, # g/cm^3
        1.04E+06   # cm
    ),
)

# Create a Monte Carlo transport engine, and configure it for backward
# transport.
engine = goupil.TransportEngine(geometry)
engine.mode = "Backward"

# Set the detector volume as a transport boundary.
engine.boundary = "Detector"

# Initialise the Monte Carlo states and source energies.
N = 1000000
states = goupil.states(N)
source_energies = numpy.empty(N)

alpha = 0.5 # This factor, which must be in [0,1], controls the fraction of
            # photo-peaks events that are simulated.
geometry.lib.initialise_states_backward(
    alpha,
    N,
    states.ctypes.data,
    source_energies.ctypes.data,
)

# Run the simulation.
status = engine.transport(states, source_energies)

# Select events consistent with a volume source located in the air.
sector = geometry.locate(states)
air_index = geometry.sector_index("Atmosphere")
selection = (status == goupil.TransportStatus.ENERGY_CONSTRAINT) & \
            (sector == air_index)
collected = states[selection]

# Print the Monte Carlo statistics.
WORLD_SIZE, DETECTOR_WIDTH, DETECTOR_HEIGHT = 2E+05, 2E+03, 1E+03
source_volume = 0.5 * WORLD_SIZE**3 - DETECTOR_WIDTH**2 * DETECTOR_HEIGHT
source_density = 1.0 / (4.0 * numpy.pi * source_volume) # A normalised source
                                                        # intensity is assumed.

rates = collected["weight"] * source_density / N
mu = sum(rates)
sigma = sum(rates**2 - (mu / N)**2)**0.5
print(f"rate = {mu:.1E} +- {sigma:.1E}")

efficiency = collected.size / N
sigma = (efficiency * (1.0 - efficiency) / N)**0.5
print(f"efficiency = {efficiency:.1E} +- {sigma:.1E}")
