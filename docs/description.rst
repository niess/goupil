Overview
========

.. _description:

----

Goupil is a Monte Carlo transport engine designed for low-energy gamma photons
(0.1-3 MeV), typically emitted by radio-isotopes like radon-222. By using a
backward algorithm and making a few approximations, Goupil can significantly
accelerate Monte Carlo simulations, resulting in computation time savings of up
to 5 orders of magnitude in certain cases. For more information, refer
to [Niess24]_.

This section provides a practical overview of Goupil illustrated with examples.
Refer to the :doc:`Python <library/index>` and :doc:`C <interface/index>`
interface sections of this document for a detailed description of the user
interfaces.

.. note::

   The following instructions assume technical familiarity with Geant4 and prior
   experience in defining a Monte Carlo geometry. If you are not familiar with
   Geant4, please refer first to the `Geant4 documentation
   <https://geant4.web.cern.ch>`_.

   In addition, it will be required to build a `shared library
   <https://en.wikipedia.org/wiki/Shared_library>`_, which is neither covered by
   this guide as there are multiple ways to achieve this depending on your OS
   and build system.


Installing goupil
-----------------

Goupil is distributed as a Python 3 module, available from PyPI. For instance,
it can be installed as:

.. code:: console

   $ pip3 install goupil

Goupil only implements transport physics. For practical applications, an
external geometry engine and software adaptation (according to the
specifications of the Goupil :doc:`C interface <interface/index>`) are also
required. Goupil comes with a source distribution of a Geant4 adapter, hereafter
called :cpp:`G4Goupil`. The corresponding source files can be accessed as:

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
:doc:`library/external_geometry` wrapper object. First, let us import
:doc:`goupil <library/index>` module as

>>> import goupil

Then, the geometry library (let's say :python:`"libgeometry.so"` on Linux) is
loaded as

>>> geometry = goupil.ExternalGeometry("path/to/libgeometry.so")

According to Goupil's model, a Monte Carlo geometry is a set of sectors that are
connected by one or more interface. Each sector is filled with a material that
has a uniform atomic composition, but its density may vary continuously.
Following, an :doc:`library/external_geometry` has two read-only attributes:
:py:attr:`materials <ExternalGeometry.materials>` and :py:attr:`sectors
<ExternalGeometry.sectors>`. These attributes list all the materials and sectors
that are defined by the geometry. For instance, as:

>>> geometry.materials
(G4_AIR, G4_SILICON_DIOXIDE)


Modifying the geometry
----------------------

The physical properties of an :doc:`library/external_geometry` can be modified
with the :py:meth:`update_material <ExternalGeometry.update_material>` and
:py:meth:`update_sector <ExternalGeometry.update_sector>` methods. For example,
let us define an exponential :doc:`library/density_gradient` to describe the air
density in the lower part of the Earth atmosphere (i.e. the troposphere).

>>> gradient = goupil.DensityGradient(1.205E-03, 1.04E+05)

Then, the density model of the first sector (index :python:`0`) can be changed
as:

>>> geometry.update_sector(0, density=gradient)

.. note::

   It is not possible to modify the structural properties of the external
   geometry, such as the number of sectors, directly from :doc:`goupil
   <library/index>`. However, it is possible to implement mutable structural
   properties at the C level in the geometry library, which can be accessed from
   Python e.g. using :external:py:mod:`ctypes`. In this case, the
   :doc:`library/external_geometry` must be reloaded whenever the Geant4
   geometry needs to be rebuilt, (i.e. when :cpp:`Geometry::Construct` is
   invoked, in the current example).


Running a simulation
--------------------

The Monte Carlo transport of photons is managed by a
:doc:`library/transport_engine` taking in charge a specific geometry. A
:doc:`library/transport_engine` is created as:

>>> engine = goupil.TransportEngine(geometry)

Each engine has its own :doc:`library/random_stream`, which can be accessed
through the :py:attr:`random <TransportEngine.random>` attribute. By default,
this stream is seeded from the system entropy. For example purposes, let us set
a specific seed value.

>>> engine.random.seed = 123456789

.. note::

   Setting a seed has the effect of reseting the pseudo-random stream.

The transport engine is set to perform a classical (forward) Monte Carlo
simulation by default. Let us instead configure the engine for backward
transport. This is done as:

>>> engine.mode = "Backward"

.. note::

   See :doc:`library/transport_settings` for a summary of configurable
   parameters.


Then, let us define a set of :python:`100` Monte Carlo states representing
photons with an energy of :python:`0.5` MeV. This is done with the
:doc:`library/states` function as

>>> states = goupil.states(100, energy=0.5)

The :doc:`library/states` function returns a `numpy structured array
<https://numpy.org/doc/stable/user/basics.rec.html>`_ of states, containing the
photons energies, their locations, etc. Since we perform a backward simulation,
these states represent final states, e.g., at a particular observation point. In
practice, one would also specify the positions and directions of observed
photons. However, for now, let us use default values for those.

Then, let us backward propagate the observed photons through the geometry. This
is done with the :py:meth:`transport <TransportEngine.transport>` method, as:

>>> status = engine.transport(states, sources_energies=1.0)

.. warning::

   The :py:meth:`transport <TransportEngine.transport>` method modifies the
   *states* array in-place. After completion, the *states* array will
   contain the propagated photons instead of the original ones.

The second argument, *sources_energies*, requires further explanation. When
running a backward Monte Carlo simulation, information about sources is needed
to correctly terminate the transport. Goupil considers two types of sources:

- Surface sources with a distributed energy spectrum, such as an external flux
  of gamma-rays.
- Volume sources with a discrete energy spectrum, such as scattered
  radio-isotopes.

In the previous example, a constant value of :python:`1.0` MeV was assumed for
the energy of volume sources.

.. note::

   The *sources_energies* argument should be omitted if there are no volume
   sources or in the case of a forward Monte Carlo.

.. note::

   In a backward transport, contained surface sources (i.e. not located on an
   outer boundary of the geometry) can be specified as a sector
   :py:attr:`boundary <TransportSettings.boundary>` at the level of the
   :doc:`library/transport_engine`.


Inspecting results
------------------

The :py:meth:`transport <TransportEngine.transport>` method returns an array of
integer codes (:doc:`library/transport_status`) which indicate the
termination condition for each propagated photon. For instance, backward
propagated photons that are consistent with a volume source can be selected as
follows:

>>> sel = (status == goupil.TransportStatus.ENERGY_CONSTRAINT)

These photons should have an energy of :python:`1.0` MeV, as requested:

>>> events[sel]["energy"]
array([1., 1., ...])

The corresponding geometry sectors can be located as:

>>> geometry.locate(events[sel])
array([1, 1, ...])

Finally, assuming a uniform volume activity of sources of :math:`1\,
\text{cm}^{-3} \text{sr}^{-1} \text{s}^{-1}`, a Monte-Carlo
estimate of the flux at the observation point is given as

>>> numpy.mean(events[sel]["weight"])

where the result has units :math:`\text{MeV}^{-1}\text{cm}^{-2} \text{sr}^{-1}
\text{s}^{-1}`.
