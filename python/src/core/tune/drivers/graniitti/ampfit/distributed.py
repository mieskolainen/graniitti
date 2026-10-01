# Persistent event partitions and exact distributed amplitude derivatives
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch
from torch.autograd.function import once_differentiable

from core.io.serialize import load_json_file
from core.tune.drivers.graniitti.ampfit.amplitude import Basis, coefficients, contract_bank


# Keep disjoint event ranges resident on a physics worker across fit evaluations
class Partition:
    # Read each assigned event range once without materializing a bank on local disk
    def __init__(self, plans, threads, batch_events):
        from core.tune.drivers.graniitti.driver import GraniittiDriver

        torch.set_num_threads(threads)
        self.banks, self.batch_events = {}, batch_events
        driver = GraniittiDriver()
        for key, directory, start, stop in plans:
            directory = Path(directory)
            metadata = load_json_file(directory / "amplitudes.bin.json")
            storage = Path(metadata.get("directory", directory))
            _, resonance, continuum = coefficients(driver, directory / "tune", load_json_file(directory / "steering.json"))
            shape = metadata["shape"]
            values = torch.empty((shape[0], stop - start, shape[2]), dtype=torch.complex128)
            mapped = (Basis(directory, metadata) if "parts" in metadata else
                      torch.from_numpy(np.memmap(storage / "amplitudes.bin", mode="c", dtype=np.complex128, shape=tuple(shape))))
            for row in range(shape[0]):
                values[row].copy_(mapped[row, start:stop])
            del mapped
            source = storage / ("source.bin" if "parts" in metadata else "amplitudes.bin")
            kinematics = torch.from_numpy(np.fromfile(str(source) + ".kinematics", dtype=np.float64,
                                                      offset=start * 5 * 8, count=(stop - start) * 5).reshape(-1, 5))
            daughters = [entry["mass"] for entry in metadata["daughters"]]
            self.banks[key] = SimpleNamespace(basis=values, resonance=resonance, continuum=continuum,
                mass2=kinematics[:, 0], transfer=kinematics[:, 1:3], daughters=daughters,
                virtuality=2 * (daughters[0]**2 - kinematics[:, 3:].amax(1)).clamp_min(0))

    # Confirm that all component reads have completed before starting optimizer trials
    def ready(self):
        return sum(bank.basis.numel() * bank.basis.element_size() for bank in self.banks.values())

    # Evaluate amplitudes, a vector Jacobian product or its exact derivative in bounded batches
    def evaluate(self, key, names, values, fixed, adjoint=None, direction=None):
        bank = self.banks[key]
        theta = torch.tensor(values, dtype=torch.float64, requires_grad=adjoint is not None)
        gradient = torch.zeros_like(theta)
        amplitudes, tangent = [], []
        for start in range(0, len(bank.mass2), self.batch_events):
            events = slice(start, start + self.batch_events)
            with torch.set_grad_enabled(adjoint is not None):
                parameters = fixed | dict(zip(names, theta, strict=True))
                amplitude = contract_bank(bank.resonance, bank.continuum, bank.basis[:, events], parameters,
                    bank.mass2[events], bank.transfer[events], bank.daughters, bank.virtuality[events])
                if adjoint is None:
                    amplitudes.append(amplitude.numpy())
                    continue
                weight = torch.tensor(adjoint[events], dtype=amplitude.dtype, requires_grad=direction is not None)
                scalar = (amplitude.conj() * weight).real.sum() + theta.sum() * 0
                first, = torch.autograd.grad(scalar, theta, create_graph=direction is not None)
                if direction is None:
                    gradient += first.detach()
                else:
                    product = (first * theta.new_tensor(direction)).sum() + theta.sum() * 0 + weight.real.sum() * 0
                    second, dual = torch.autograd.grad(product, (theta, weight))
                    gradient += second.detach()
                    tangent.append(dual.numpy())
        if adjoint is None:
            return np.concatenate(amplitudes)
        return gradient.numpy() if direction is None else (gradient.numpy(), np.concatenate(tangent))


# Preserve exact first and second derivatives while exchanging only contracted event amplitudes
class RemoteAmplitude(torch.autograd.Function):
    # Evaluate independent cached partitions on their assigned physics workers
    @staticmethod
    def forward(ctx, theta, workers, names, fixed):
        import ray

        ctx.save_for_backward(theta)
        ctx.workers, ctx.names, ctx.fixed = workers, names, fixed
        arrays = ray.get([actor.evaluate.remote(key, names, theta.detach().numpy(), fixed)
                          for actor, key, _, _ in workers])
        return torch.from_numpy(np.concatenate(arrays))

    # Recompute local graphs from resident components without retaining a full bank graph
    @staticmethod
    def backward(ctx, adjoint):
        theta, = ctx.saved_tensors
        return RemoteVJP.apply(theta, adjoint, ctx.workers, ctx.names, ctx.fixed), None, None, None


# Differentiate the distributed vector Jacobian product for covariance calculations
class RemoteVJP(torch.autograd.Function):
    # Sum worker gradients in the physical event partition order
    @staticmethod
    def forward(ctx, theta, adjoint, workers, names, fixed):
        import ray

        ctx.save_for_backward(theta, adjoint)
        ctx.workers, ctx.names, ctx.fixed = workers, names, fixed
        arrays = ray.get([actor.evaluate.remote(key, names, theta.detach().numpy(), fixed,
                                               adjoint[start:stop].detach().resolve_conj().resolve_neg().numpy())
                          for actor, key, start, stop in workers])
        return torch.from_numpy(np.sum(arrays, axis=0))

    # Compute both Hessian and adjoint derivatives on cached partitions
    @staticmethod
    @once_differentiable
    def backward(ctx, direction):
        import ray

        theta, adjoint = ctx.saved_tensors
        pairs = ray.get([actor.evaluate.remote(key, ctx.names, theta.detach().numpy(), ctx.fixed,
                            adjoint[start:stop].detach().resolve_conj().resolve_neg().numpy(), direction.detach().numpy())
                         for actor, key, start, stop in ctx.workers])
        return (torch.from_numpy(np.sum([pair[0] for pair in pairs], axis=0)),
                torch.from_numpy(np.concatenate([pair[1] for pair in pairs])), None, None, None)


# Bind differentiable coordinates independently of fixed steering values
def remote_amplitude(parameters, workers):
    names = sorted(name for name, value in parameters.items() if torch.is_tensor(value))
    theta = torch.stack([parameters[name] for name in names]) if names else torch.empty(0, dtype=torch.float64)
    fixed = {name: value for name, value in parameters.items() if name not in names}
    return RemoteAmplitude.apply(theta, workers, names, fixed)


# Partition the total component bytes across admitted worker nodes before scheduling trials
def start_workers(driver, param):
    import ray
    from ray.util.scheduling_strategies import NodeAffinitySchedulingStrategy

    nodes = sorted((node for node in ray.nodes() if node.get("Alive") and
                    node["Resources"].get("CPU", 0) >= 2 and not node["Resources"].get("icetune_head", 0)),
                   key=lambda node: node["NodeID"])
    if not nodes:
        raise ValueError("Distributed amplitude banks require worker nodes with at least two CPUs")
    directory = driver.amplitude_directory(cdir=param["cdir"], run_name=param["run_name"])
    samples = []
    for path in sorted(directory.glob("*/*/amplitudes.bin.json")):
        key = path.parent.relative_to(directory).as_posix()
        metadata = load_json_file(path)
        storage = str(Path(metadata.get("directory", path.parent)).resolve())
        columns, events, helicities = metadata["shape"]
        samples.append((key, storage, events, columns * helicities * np.dtype(np.complex128).itemsize))
    if not samples:
        raise ValueError("No prepared amplitude samples were found")
    total = sum(events * size for _, _, events, size in samples)
    capacity = sum(node["Resources"].get("memory", 0) for node in nodes)
    if total * 3 > capacity:
        raise ValueError(f"Amplitude partitions require {total * 3} bytes of worker memory, available {int(capacity)}")
    plans = [[] for _ in nodes]
    sample_index, start = 0, 0
    for index, node in enumerate(nodes):
        budget = math.ceil(total * node["Resources"]["memory"] / capacity)
        used = 0
        while sample_index < len(samples) and (used < budget or index == len(nodes) - 1):
            key, storage, events, size = samples[sample_index]
            count = min(events - start, max(1, (budget - used) // size)) if index < len(nodes) - 1 else events - start
            plans[index].append((key, storage, start, start + count))
            used += count * size
            start += count
            if start == events:
                sample_index += 1
                start = 0
    workers, ready = {}, []
    actors = []
    try:
        for node, plan in zip(nodes, plans, strict=True):
            if not plan:
                continue
            size = sum((stop - start) * next(sample[3] for sample in samples if sample[0] == key)
                       for key, _, start, stop in plan)
            # Cache, loading pages and differentiation buffers each fit within one partition size
            memory = size * 3
            if memory > node["Resources"]["memory"]:
                raise ValueError("An event partition exceeds the admitted worker memory")
            cpus = int(node["Resources"]["CPU"]) - 1
            actor = ray.remote(Partition).options(num_cpus=cpus, memory=memory, max_restarts=-1,
                max_task_retries=-1, scheduling_strategy=NodeAffinitySchedulingStrategy(node["NodeID"], soft=True)).remote(
                    plan, cpus, param["mc_steer"]["ampfit"]["batch_events"])
            actors.append(actor)
            ready.append(actor.ready.remote())
            for key, _, start, stop in plan:
                workers.setdefault(key, []).append((actor, key, start, stop))
        print(f"ampfit cache: loading {total / 1024**3:.4f} GiB on {len(actors)} worker nodes", flush=True)
        ray.get(ready, timeout=param["max_t"])
        print("ampfit cache: all event partitions are resident", flush=True)
    except BaseException:
        for actor in actors:
            ray.kill(actor)
        raise
    return workers

