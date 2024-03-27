#! /usr/bin/env python3
import ctypes
import goupil


# Load the external geometry (using an explicit relative path).
geometry = goupil.ExternalGeometry("./libgeometry.so")

# Configure the user defined external function(s) (using ctypes).
geometry.lib.initialise_states_forward.argtypes = [
    ctypes.c_size_t,
    ctypes.c_void_p,
]
geometry.lib.initialise_states_forward.restype = None

# Set a vertical density gradient for the atmosphere.
geometry.update_sector(
    sector = "Atmosphere",
    density = goupil.DensityGradient(
        1.225E-03, # g/cm^3
        1.04E+06   # cm
    )
)

# Create a Monte Carlo transport engine.
engine = goupil.TransportEngine(geometry)

# Set the detector volume as a transport boundary.
engine.boundary = "Detector"

# Initialise the Monte Carlo states.
N = 1000000
states = goupil.states(N)
geometry.lib.initialise_states_forward(N, states.ctypes.data)

# Run the simulation.
status = engine.transport(states)

# Select collected events.
collected = states[status == goupil.TransportStatus.BOUNDARY]

# Print the Monte Carlo statistics. Note that in this case the collection rate
# and Monte Carlo efficiency are proportional.

source_density = 1E-05 # Bq/cm^3
WORLD_SIZE, DETECTOR_SIZE = 2E+05, 2E+03 # cm
source_volume = 0.5 * (WORLD_SIZE**3 - DETECTOR_SIZE**3)
total_activity = source_density * source_volume

efficiency = collected.size / N
sigma_efficiency = (efficiency * (1.0 - efficiency) / N)**0.5

rate = efficiency * total_activity * 1E-06 # MHz
sigma_rate = sigma_efficiency * total_activity * 1E-06 # MHz

print(f"rate = {rate:.2E} +- {sigma_rate:.2E} MHz")
print(f"efficiency = {efficiency:.1E} +- {sigma_efficiency:.1E}")
