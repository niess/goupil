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
