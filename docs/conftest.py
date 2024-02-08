import goupil
import numpy
import pytest

numpy.set_printoptions(precision=5)

@pytest.fixture(autouse=True)
def add_goupil(doctest_namespace):
    doctest_namespace["goupil"] = goupil
    doctest_namespace["numpy"] = numpy
