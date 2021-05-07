# Pytest configuration file

import pytest

def pytest_addoption(parser):
    parser.addoption("--POMLOOP", type=str, default="false")
    parser.addoption("--NEVENTS", type=int, default=50000)
    parser.addoption("--MODELPARAM", type=str, default="TUNE0")

@pytest.fixture(scope='session')
def POMLOOP(request):
    value = request.config.option.POMLOOP
    value = True if value.lower() == "true" else False
    
    if value is None:
        pytest.skip()
    return value

@pytest.fixture(scope='session')
def NEVENTS(request):
    value = request.config.option.NEVENTS
    if value is None:
        pytest.skip()
    return value

@pytest.fixture(scope='session')
def MODELPARAM(request):
    value = request.config.option.MODELPARAM
    if value is None:
        pytest.skip()
    return value
