Overview
========

.. _description:

----

Goupil is a Monte Carlo transport engine designed for low-energy gamma photons
(0.1-3 MeV), typically emitted by radio-isotopes like radon-222. By using a
backward algorithm and making a few approximations, Goupil significantly
accelerates Monte Carlo transport simulation, resulting in computation time
savings of up to 5 orders of magnitude in certain cases. For more information,
refer to [Niess24]_.

This section provides a practical overview of Goupil and includes examples of
its use. Refer to the :doc:`Python <library/index>` and :doc:`C
<interface/index>` interface sections of this document for a detailed
description of the user interfaces.

.. note::

   The following instructions assume technical familiarity with Geant4 and prior
   experience in defining a Monte Carlo geometry. If you are not familiar with
   Geant4, please refer to the `Geant4 documentation
   <https://geant4.web.cern.ch>`_ for instructions on installation and defining
   a geometry.

   In addition, it will be required to build a `shared library
   <https://en.wikipedia.org/wiki/Shared_library>`_, which is neither covered by
   this guide as there are multiple ways to achieve this depending on your OS
   and build system.


Installation
------------

Goupil is distributed as a Python 3 module, available from PyPI. For instance,
it can be installed as follows.

.. code:: console

   $ pip3 install goupil

Goupil only implements transport physics. For practical applications, an
external geometry engine and software adaptation (according to the
specifications of the Goupil :doc:`C interface <interface/index>`) are also
required. Goupil comes with a source distribution of an adapter to Geant4,
hereafter called :cpp:`G4Goupil`. These source files can be accessed as follows.

.. code:: console

   $ ls "$(python3 -m goupil --prefix)/interfaces/geant4"
   G4Goupil.cc  G4Goupil.hh  goupil.h

.. warning::

   Using :cpp:`G4Goupil` requires a working installation of Geant4.


Preparing the geometry
----------------------

:cpp:`G4Goupil` allows a Geant4 geometry to be exported as a shared library,
which can then be navigated using Goupil's Python interface. To illustrate this
mechanism, consider the following example. Let :cpp:`Geometry` be a subclass of
:cpp:`G4VDetectorConstruction` that implements a Geant4 geometry. Then, embed
the following :cpp:`G4Goupil::NewGeometry` function in the geometry library.

.. code:: C++

   const G4VPhysicalVolume * G4Goupil::NewGeometry() {
       // Build the geometry and return the top "World" volume.
       Geometry geometry();
       return geometry.Construct();
   }

.. warning::

   The geometry shared library must link to or include both the :cpp:`Geometry`
   and :cpp:`G4Goupil` implementations. This means that G4Goupil.cc must be
   compiled at some point.

Optionally, a cleanup function (:cpp:`G4Goupil::DropGeometry`) can be defined
for when the geometry is released by Goupil. The implementation will depend on
the memory policy used when building the geometry. For example, if materials are
managed by a global store (e.g. :cpp:`G4NistManager`), the following code could
be used.

.. code:: C++

   void G4Goupil::DropGeometry(const G4VPhysicalVolume * volume) {
       // Delete any sub-volume(s).
       auto && logical = volume->GetLogicalVolume();
       while (logical->GetNoDaughters()) {
           auto daughter = logical->GetDaughter(0);
           logical->RemoveDaughter(daughter);
           G4Goupil::DropGeometry(daughter);
       }
       // Delete this volume.
       delete logical->GetSolid();
       delete logical;
       delete volume;
   }


Importing the geometry
----------------------

The geometry library is imported using an :doc:`library/external_geometry`
wrapper object. First, let us import :doc:`goupil <library/index>` module as

>>> import goupil

Then, the geometry library (let's say :python:`"libgeometry.so"` on Linux) is
loaded as

>>> geometry = goupil.ExternalGeometry("path/to/libgeometry.so")

According to Goupil, a Monte Carlo geometry is a set of sectors that are
connected by interfaces. Each sector is filled with a material that has a
uniform atomic composition, but its density may vary continuously. The
:doc:`library/external_geometry` object has two read-only attributes:
:py:attr:`materials <ExternalGeometry.materials>` and :py:attr:`sectors
<ExternalGeometry.sectors>`. These attributes list all the materials and sectors
that are defined by the loaded geometry. For instance,

>>> geometry.materials
(G4_AIR, G4_SILICON_DIOXIDE)


Modifying the geometry
----------------------

To modify the physical properties of geometry sectors, the
:doc:`library/external_geometry` object provides the :py:meth:`update_material
<ExternalGeometry.update_material>` and :py:meth:`update_sector
<ExternalGeometry.update_sector>` methods. For example, let us define an
exponential :doc:`library/density_gradient` to describe the air density in the
lower atmosphere (i.e. the troposphere).

>>> gradient = goupil.DensityGradient(1.205E-03, 1.04E+05)

Then, the density model of the first sector (index :python:`0`) can be changed
as follows.

>>> geometry.update_sector(0, density=gradient)

.. note::

   It is not possible to modify the structural properties of the external
   geometry, such as the number of sectors, from :doc:`goupil <library/index>`.
   However, it is possible to implement mutable structural properties at the C
   level in the geometry library, which can be accessed from Python e.g. using
   :external:py:mod:`ctypes`. In this case, the :doc:`library/external_geometry`
   must be reloaded whenever the Geant4 geometry needs to be rebuilt, (i.e. when
   :cpp:`Geometry::Construct` is invoked, in the current example).
