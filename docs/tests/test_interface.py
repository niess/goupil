import goupil
import numpy
import pytest


def test_AtomicElement():
    """Test usage of an AtomicElement."""

    # Check constructor from Z.
    for z in range(1, 101):
        element = goupil.AtomicElement(z)
        assert element.Z == z
        assert (element.A > 0.0) and (element.A != z)

    # Check out of range values.
    for z in (0, 119):
        with pytest.raises(RuntimeError) as e:
            element = goupil.AtomicElement(z)
        assert str(e.value).startswith("bad atomic number")

    # Check constructor from symbol.
    for symbol in ("H", "C", "Fe", "U"):
        element = goupil.AtomicElement(symbol)
        assert element.symbol == symbol
        assert isinstance(element.name, str)
        assert isinstance(element.A, float)
        assert isinstance(element.Z, int)

    # Check unknown symbol.
    with pytest.raises(RuntimeError) as e:
        element = goupil.AtomicElement("Zx")
    assert str(e.value).startswith("no such atomic element")

    # Check comparison.
    H0 = goupil.AtomicElement("H")
    H1 = goupil.AtomicElement("H")
    assert H0 == H1

    # Check immutable.
    for attr in ("A", "name", "symbol", "Z"):
        with pytest.raises(AttributeError) as e:
            setattr(H0, attr, None)
        assert "not writable" in str(e.value)


def test_ComptonProcess():
    """Test usage of a ComptonProcess."""

    # Check constructor.
    process = goupil.ComptonProcess()
    assert process.method == "Rejection Sampling"
    assert process.mode == "Direct"
    assert process.model == "Scattering Function"
    assert process.precision == 1.0

    for method in ("Inverse Transform", "Rejection Sampling"):
        for mode in ("Adjoint", "Direct", "Inverse"):
            for model in ("Klein-Nishina", "Penelope", "Scattering Function"):
                try:
                    process = goupil.ComptonProcess(
                        method=method,
                        mode=mode,
                        model=model
                    )
                except NotImplementedError as e:
                    if str(e).startswith("bad sampling"):
                        continue
                else:
                    assert process.method == method
                    assert process.mode == mode
                    assert process.model == model
                    assert process.precision == 1.0

    process = goupil.ComptonProcess(precision=10.0)
    assert process.precision == 10.0

    with pytest.raises(ValueError):
        goupil.ComptonProcess(precision=0)

    with pytest.raises(KeyError):
        goupil.ComptonProcess(toto=0)

    # Check cross-section method.
    process = goupil.ComptonProcess()
    H = goupil.elements("H")
    material = goupil.MaterialDefinition(
        name = "Material",
        mole_composition = ((1, H),)
    )
    assert process.cross_section(1.0, material) > 0.0

    process = goupil.ComptonProcess(model="Klein-Nishina")
    energies = numpy.logspace(-2, 1, 31)
    values = process.cross_section(energies, material)
    assert values.shape == energies.shape
    assert (numpy.diff(values) < 0.0).all()


def test_MaterialDefinition():
    """Test usage of a MaterialDefinition."""

    # Check constructor.
    H2O = goupil.MaterialDefinition("H2O")
    composition = H2O.mole_composition
    assert len(composition) == 2
    assert composition[0][0] == 2
    assert str(composition[0][1]) == "H"
    assert composition[1][0] == 1
    assert str(composition[1][1]) == "O"
    assert H2O.name == "H2O"

    H, O = goupil.elements("H, O")
    assert H2O.mass == 2 * H.A + O.A

    nothing = goupil.MaterialDefinition()
    assert nothing.mass == 0.0

    water = goupil.MaterialDefinition(
        name = "Water",
        mole_composition = (
            (2, "H"),
            (1, O)
        )
    )
    assert water.name == "Water"
    assert water.mass == H2O.mass
    assert water.mole_composition == H2O.mole_composition

    mixture = goupil.MaterialDefinition(
        name = "Mixture",
        mass_composition = (
            (0.5, H2O),
            (0.5, water)
        )
    )
    assert mixture.mass == water.mass
    assert mixture.mole_composition == water.mole_composition

    with pytest.raises(RuntimeError) as e:
        goupil.MaterialDefinition("Xu")
    assert str(e.value) == "no such atomic element 'Xu'"


def test_MaterialRecord():
    """Test usage of a MaterialRecord."""

    # Check direct instanciation.
    with pytest.raises(TypeError) as e:
        goupil.MaterialRecord()
    assert str(e.value) == "No constructor defined"

    # Check attributes.
    H2O = goupil.MaterialDefinition("H2O")
    registry = goupil.MaterialRegistry(H2O)
    registry.compute()

    record = registry["H2O"]
    assert(record.definition == H2O)
    assert(record.electrons == H2O.electrons())

    # Check table getters.
    table = record.absorption_cross_section()
    assert(isinstance(table, goupil.CrossSection))
