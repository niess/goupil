def pytest_collection_modifyitems(session, config, items):
    """Monkey patch tests' order."""

    def rank(a):
        return \
            2**2 * int("requires_matplotlib" in a.keywords) + \
            2**1 * int("requires_calzone" in a.keywords) + \
            2**0 * int("example" in a.keywords)

    items.sort(key=rank)
