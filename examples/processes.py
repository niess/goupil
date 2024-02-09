#! /usr/bin/env python3
"""
This example demonstrates how to examine physical properties of the interaction
processes implemented in Goupil. The following modules will be used for this
purpose.
"""

import goupil
import numpy
import matplotlib.pyplot as plot


"""
Goupil implements three types of interaction processes between gamma rays and
matter: absorption, Compton scattering, and Rayleigh scattering. The absorption
process includes the photoelectric effect and e+e- pair production. Let us
define interfaces for these processes as follows.
"""

absorption = goupil.AbsorptionProcess
compton = goupil.ComptonProcess()
rayleigh = goupil.RayleighProcess

"""
We will start by comparing the total cross-sections of these processes,
using water as the target material.
"""

energies = numpy.logspace(-3, 1, 161)
material = "H2O"

plot.figure()
plot.loglog(energies, absorption.cross_section(energies, material), "k--",
            label="Absorption")
plot.loglog(energies, compton.cross_section(energies, material), "k-",
            label="Compton")
plot.loglog(energies, rayleigh.cross_section(energies, material), "k:",
            label="Rayleigh")
plot.xlabel("energy (MeV)")
plot.ylabel("cross-section (cm$^2$)")
plot.legend()

"""
The Compton process can be modelled in several ways. A simplified approach is to
consider atomic electrons as free and at rest, which leads to the Klein-Nishina
cross section. As an example, we will compare this model to the one used by
default by Goupil (i.e. the Scattering Function model).
"""

klein_nishina = goupil.ComptonProcess(model="Klein-Nishina")

plot.figure()
plot.semilogx(energies, klein_nishina.cross_section(energies, material), "k--",
              label="Klein-Nishina")
plot.semilogx(energies, compton.cross_section(energies, material), "k-",
              label="Scattering Function")
plot.xlabel("energy (MeV)")
plot.ylabel("cross-section (cm$^2$)")
plot.legend()

"""
The Penelope model is also included for comparison purposes. However, it is not
applicable in backward Monte Carlo simulations. Goupil's default model, by
design, yields the same total effective cross-section as Penelope's. We will now
compare the differential cross sections of these three models, while considering
a gold target.
"""

material = "Au"
incoming_energy = 0.05 # MeV
outgoing_energies = numpy.linspace(0.7 * incoming_energy, incoming_energy, 101)

plot.figure()
plot.plot(outgoing_energies / incoming_energy,
          klein_nishina.dcs(incoming_energy, outgoing_energies, material),
          "k--", label="Klein-Nishina")
plot.plot(outgoing_energies / incoming_energy,
          compton.dcs(incoming_energy, outgoing_energies, material),
          "k-", label="Scattering Function")

"""
To determine the differential cross-section with respect to energy, the Penelope
model necessitates numerical integration. The Monte Carlo technique will be used
to carry out this integration.
"""

penelope = goupil.ComptonProcess(model="Penelope")
sampled_energies, *_ = penelope.sample(numpy.full(10000000, incoming_energy),
                                       material)
penelope_pdf, edges = numpy.histogram(sampled_energies, outgoing_energies,
                                      density=True)

x = 0.5 * (edges[1:] + edges[:-1]) / incoming_energy
y = penelope_pdf * penelope.cross_section(incoming_energy, material)

plot.plot(x, y, "r-", label="Penelope")
plot.xlabel("energy (MeV)")
plot.ylabel("dcs (cm$^2$/MeV)")
plot.legend()

plot.show()
