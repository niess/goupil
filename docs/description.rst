Description
===========

.. _description:

----

Goupil is a Monte Carlo transport engine designed for low-energy gamma photons
(0.1-3 MeV), typically emitted by radio-isotopes like radon-222. By using a
backward algorithm and making a few approximations, Goupil significantly
accelerates Monte Carlo transport simulation, resulting in computation time
savings of up to 5 orders of magnitude in certain cases. For more information,
refer to [Niess24]_.

This section provides a practical description of Goupil and includes examples of
its use. Refer to the :doc:`Python <library/index>` and :doc:`C
<interface/index>` interface sections of this document for a detailed
description of the user interfaces.


Installation
------------

Goupil is distributed as a Python 3 module, available from PyPI. For instance,
it can be installed as follows:

.. code:: console

   $ pip3 install goupil

Goupil implements only transport physics. For practical applications, an
external geometry engine and software adaptation (according to the
specifications of the Goupil :doc:`C interface <interface/index>`) are also
necessary. Goupil comes with an adapter to Geant4, accessible (after
installation of the Python module) as

.. code:: console

   $ ls $(goupil-config --prefix)/interfaces/geant4
   G4Goupil.cc  G4Goupil.hh

In addition, it is also necessary to have a working installation of Geant4 in
this case. If required, refer to the `official documentation
<https://geant4.web.cern.ch>`_ for instructions on how to install Geant4.


Examples
--------

Geant4 geometry
~~~~~~~~~~~~~~~

.. code:: C++

   const G4VPhysicalVolume * G4Goupil::NewGeometry() {
       // Build the geometry and return the top "World" volume.
       DetectorConstruction geometry();
       return geometry.Construct();
   }

.. code:: C++

   void G4Goupil::DropGeometry(const G4VPhysicalVolume * volume) {
       // Delete any sub-volumes.
       auto logical = volume->GetLogicalVolume();
       const int n = logical->GetNoDaughters();
       for (int i = 0; i < n; i++) {
           G4Goupil::DropGeometry(logical->GetDaughter(i));
       }
       // Delete this volume.
       delete logical->GetSolid();
       delete logical;
       delete volume;
   }
