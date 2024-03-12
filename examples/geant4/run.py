#! /usr/bin/env python3
import ctypes
import goupil


# Load the external geometry (using an explicit relative path).
geometry = goupil.ExternalGeometry("./libgeometry.so")

# Configure the user defined external function(s) (using ctypes).
geometry.lib.initialise_states.argtypes = [
    ctypes.c_size_t,
    ctypes.c_void_p,
    ctypes.c_bool
]
geometry.lib.initialise_states.restype = None

# Create a Monte Carlo transport engine.
engine = goupil.TransportEngine(geometry)

# Locate the detector sector, and set it as a transport boundary.
sector_names = [sector.description for sector in geometry.sectors]
detector_index = sector_names.index("Detector")
engine.boundary = detector_index

# Initialise the Monte Carlo states.
states = goupil.states(1000000)
geometry.lib.initialise_states(states.size, states.ctypes.data, True)

# Run the simulation.
status = engine.transport(states)

# Select collected events.
collected = states[status == goupil.TransportStatus.BOUNDARY]

# Print the Monte Carlo efficiency.
efficiency = collected.size / states.size
sigma = collected.size**0.5 / states.size
print(f"efficiency = {efficiency:.1E} +- {sigma:.1E}")
