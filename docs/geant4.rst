Geant4 adapter
==============

.. _geant4:

This section provides specific instructions for importing a Geant4 geometry in
Goupil. Refer to the :doc:`overview` section of this document for a general
description of Goupil.

.. note::

   The following instructions assume technical familiarity with Geant4 and prior
   experience in defining a Monte Carlo geometry. If you are not familiar with
   Geant4, please refer first to the `Geant4 documentation
   <https://geant4.web.cern.ch>`_.

   In addition, it will be required to build a `shared library
   <https://en.wikipedia.org/wiki/Shared_library>`_, which is neither covered by
   this guide as there are multiple ways to achieve this depending on your OS
   and build system.


Locating G4Goupil
-----------------

Goupil comes with a source distribution of a Geant4 adapter, hereafter called
:cpp:`G4Goupil`. The corresponding source files can be accessed as:

.. code:: console

   $ ls "$(python3 -m goupil --prefix)/interfaces/geant4"
   G4Goupil.cc  G4Goupil.hh  goupil.h

.. warning::

   Using :cpp:`G4Goupil` requires a working installation of Geant4.


Preparing the geometry
----------------------

:cpp:`G4Goupil` allows a Geant4 geometry to be exported as a shared library,
which can then be navigated using Goupil's Python interface. To illustrate this
mechanism, let us consider the following example. Let :cpp:`Geometry` be a
subclass of :cpp:`G4VDetectorConstruction` that implements a Geant4 geometry.
This geometry is exported by embedding the following
:cpp:`G4Goupil::NewGeometry` function in a shared library.

.. code:: C++

   const G4VPhysicalVolume * G4Goupil::NewGeometry() {
       // Build the geometry and return the top "World" volume.
       Geometry geometry();
       return geometry.Construct();
   }

.. warning::

   The geometry shared library must also link to or include both the
   :cpp:`Geometry` and :cpp:`G4Goupil` implementations. In particular, this
   means that G4Goupil.cc must be compiled at some point.

Optionally, a cleanup function (:cpp:`G4Goupil::DropGeometry`) can also be
included in the shared library, for when the geometry is released by Goupil. The
implementation of this function must be consistent with the memory policy used
when building the geometry. For example, assuming that materials are managed by
a global store (e.g. :cpp:`G4NistManager`), the following code could be used.

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

The previous geometry library is imported in Python by using an
:doc:`py/external_geometry` wrapper object. First, let us import :doc:`goupil
<py/index>` module as

>>> import goupil

Then, the geometry library (let's say :python:`"libgeometry.so"` on Linux) is
loaded as

>>> geometry = goupil.ExternalGeometry("path/to/libgeometry.so") # doctest: +SKIP

According to Goupil's model, a Monte Carlo geometry is a set of sectors that are
connected by one or more interface. Each sector is filled with a material that
has a uniform atomic composition, but its density may vary continuously.
Following, an :doc:`py/external_geometry` has two read-only attributes:
:py:attr:`materials <ExternalGeometry.materials>` and :py:attr:`sectors
<ExternalGeometry.sectors>`. These attributes list all the materials and sectors
that are defined by the geometry. For instance, as:

>>> geometry.materials # doctest: +SKIP
(G4_AIR, G4_CALCIUM_CARBONATE)


Modifying the geometry
----------------------

The physical properties of an :doc:`py/external_geometry` can be modified with
the :py:meth:`update_material <ExternalGeometry.update_material>` and
:py:meth:`update_sector <ExternalGeometry.update_sector>` methods. For example,
let us define an exponential :doc:`py/density_gradient` to describe the air
density in the lower part of the Earth atmosphere (i.e. the troposphere).

>>> gradient = goupil.DensityGradient(1.205E-03, 1.04E+05) # doctest: +SKIP

Then, the density model of the first sector (index :python:`0`) can be changed
as:

>>> geometry.update_sector(0, density=gradient) # doctest: +SKIP

.. note::

   It is not possible to modify the structural properties of the external
   geometry, such as the number of sectors, directly from :doc:`goupil
   <py/index>`. However, it is possible to implement mutable structural
   properties at the C level in the geometry library, which can be accessed from
   Python e.g. using :external:py:mod:`ctypes`. In this case, the
   :doc:`py/external_geometry` must be reloaded whenever the Geant4 geometry
   needs to be rebuilt, (i.e. when :cpp:`Geometry::Construct` is invoked, in the
   current example).
