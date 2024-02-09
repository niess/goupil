#! /usr/bin/env python3
import goupil
import numpy
import matplotlib.pyplot as plot


energies = numpy.logspace(-3, 1, 161)
material = "H2O"
absorption = goupil.AbsorptionProcess.cross_section(energies, material)
compton = goupil.ComptonProcess().cross_section(energies, material)
rayleigh = goupil.RayleighProcess.cross_section(energies, material)

plot.figure()
plot.loglog(energies, absorption, "k--", label="Absorption")
plot.loglog(energies, compton, "k-", label="Compton")
plot.loglog(energies, rayleigh, "k:", label="Rayleigh")
plot.xlabel("energy (MeV)")
plot.ylabel("cross-section (cm$^2$)")
plot.legend()
plot.show()
