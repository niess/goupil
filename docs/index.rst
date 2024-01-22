Goupil 0.2.0
============

Goupil is a Monte Carlo transport engine designed for low-energy gamma photons
(0.1-3 MeV), typically emitted by radio-isotopes like radon-222 descendants. By
using a backward algorithm and making a few approximations, Goupil can
significantly accelerate Monte Carlo simulations, resulting in computation time
savings of up to 5 orders of magnitude in certain cases. For more information,
or for citing Goupil, please refer to [Niess24]_.


Documentation
-------------

.. toctree::
   :maxdepth: 1

   overview
   py/index
   c/index
   references


System of units
---------------

.. note::

   Goupil uses the Centimetre-Gram-Second (CGS) system of units (e.g.
   g/cm\ :sup:`3` for a density), except for energies and momenta
   which are expressed in MeV and MeV/c respectively.
