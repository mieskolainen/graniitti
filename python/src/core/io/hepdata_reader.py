# Load HEPData readers declared by icepack bundles
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import inspect

from core.io import steering


# Load the reader declared by one icepack bundle
def load_reader(reference: str, *, dataset_path: str, cdir: str):
    resolved = steering.resolve_python_reference(
        reference, package="core.io.hepdata_reader", dataset_path=dataset_path, cdir=cdir
    )
    module = steering.load_python_module(resolved)
    if not callable(getattr(module, "read", None)):
        raise AttributeError(f"HEPData reader '{resolved}' must define callable read")
    return module


# Load the reader declared inside one dataset card
def load_dataset_reader(dataset_reference: str, *, cdir: str):
    dataset, dataset_path = steering.load_dataset(dataset_reference, cdir=cdir)
    return load_reader(dataset["reader"], dataset_path=dataset_path, cdir=cdir)


# Read a table through the reader declared by its icepack bundle
def read(reference: str, *, dataset_path: str, cdir: str, **kwargs) -> dict:
    reader = load_reader(reference, dataset_path=dataset_path, cdir=cdir).read
    signature = inspect.signature(reader)
    if any(parameter.kind == inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values()):
        return reader(**kwargs)
    accepted = {key: value for key, value in kwargs.items() if key in signature.parameters}
    return reader(**accepted)
