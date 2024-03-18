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

# Create a Monte Carlo transport engine.
engine = goupil.TransportEngine(geometry)

# Locate the detector sector, and set it as a transport boundary.
sector_names = [sector.description for sector in geometry.sectors]
detector_index = sector_names.index("Detector")
engine.boundary = detector_index

# Initialise the Monte Carlo states.
states = goupil.states(1000000)
geometry.lib.initialise_states_forward(states.size, states.ctypes.data)

# Run the simulation.
status = engine.transport(states)

# Select collected events.
collected = states[status == goupil.TransportStatus.BOUNDARY]

# Print the Monte Carlo statistics. Note that in this case the Monte Carlo
# efficiency equals the normalised rate of collected events.
m, n = collected.size, states.size
efficiency = m / n
sigma = ((1.0 - m / n) / n)**0.5
print(f"efficiency / rate = {efficiency:.1E} +- {sigma:.1E}")
