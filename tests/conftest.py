# Pytest configuration file

import pytest

def pytest_addoption(parser):
    parser.addoption("--POMLOOP", action="store", default=True)
    parser.addoption("--NEVENTS", action="store", default=50000)

@pytest.fixture(scope='session')
def POMLOOP(request):
    value = request.config.option.POMLOOP
    if value is None:
        pytest.skip()
    return value

@pytest.fixture(scope='session')
def NEVENTS(request):
    value = request.config.option.NEVENTS
    if value is None:
        pytest.skip()
    return value
