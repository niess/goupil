import goupil
import pytest

@pytest.fixture(autouse=True)
def add_goupil(doctest_namespace):
    doctest_namespace["goupil"] = goupil
