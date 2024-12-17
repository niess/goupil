Goupil
======

(**G**\ amma transp\ **O**\ rt **U**\ tility, a\ **P**\ proximate but revers\ **I**\ b\ **L**\ e)

Goupil is a Monte Carlo transport engine designed for low-energy gamma-rays
(0.1-3 MeV), typically emitted from radioactive isotopes. By using a backward
transport algorithm, and making a few approximations, Goupil can significantly
accelerate Monte Carlo simulations, resulting in computation time savings of up
to 5 orders of magnitude in the air. For an academic description of Goupil, or
for citing Goupil, please refer to [NVT24]_.


Documentation
-------------

.. toctree::
   :maxdepth: 1

   overview
   geant4
   py/index
   c/index
   references
   Source code <https://github.com/niess/goupil>


System of units
---------------

.. note::

   Goupil uses the Centimetre-Gram-Second (CGS) system of units (e.g.
   g/cm\ :sup:`3` for a density), except for energies and momenta
   which are expressed in MeV and MeV/c respectively.
