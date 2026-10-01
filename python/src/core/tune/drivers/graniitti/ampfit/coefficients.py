# Differentiable resonance coefficients from GRANIITTI tuning coordinates for ampfit
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path

import pyjson5
import torch

from core.io import jsonref
from core.io.serialize import load_json_file
from core.numerics import array
from core.numerics.interp import ChebyshevGrid
from core.tune.drivers.graniitti.tunesetup.continuum import _continuum_entry_for_pair, _production_entry_for_pair
from core.tune.drivers.graniitti.tunesetup.domains import angular_row_name, mp_spin_sectors


# Compute a unit complex phase while preserving the selected array backend
def phase(value):
    value = array.asarray(value, dtype=float)
    xp = array.namespace(value)
    return xp.cos(value) + 1j * xp.sin(value)


# Interpolate the screened residual after extracting the external power transfer factors
def transfer_weights(axis, value, transfer, power, inverse):
    nodes = axis.nodes.to(transfer.device)
    value = torch.as_tensor(value, dtype=nodes.dtype, device=transfer.device)
    cardinal = axis(value)
    if not inverse:
        nodes, value = nodes.reciprocal(), value.reciprocal()
    ratio = ((1 - nodes[None, :, None] * transfer[:, None, :]) / (1 - value * transfer[:, None, :]))
    return cardinal * ratio.pow(power).prod(-1)


# Compute the two body normalization with the native width steps below threshold
def decay_phase_space(mass, width, daughters, controls):
    mass = array.asarray(mass, like=width, dtype=float)
    width = array.asarray(width, like=mass)
    xp = array.namespace(mass)
    first, second = daughters
    effective = mass
    found = xp.zeros_like(mass) > 0.0
    for step in range(controls["decay_steps"] + 1):
        candidate = mass + step * width
        beta2 = (1 - ((first + second) / candidate)**2) * (1 - ((first - second) / candidate)**2)
        # Clamp only the branch selection, leaving the chosen physical phase space differentiable
        trial = array.sqrt(xp.maximum(beta2, array.asarray(0.0, like=mass))) / (16 * math.pi * candidate)
        valid = (candidate > first + second) & (trial > controls["decay_zero"])
        effective = xp.where(~found & valid, candidate, effective)
        found = found | valid
    beta2 = (1 - ((first + second) / effective)**2) * (1 - ((first - second) / effective)**2)
    return xp.where(found, array.sqrt(beta2) / (16 * math.pi * effective), array.asarray(float("nan"), like=mass))


# Decode the same physical parameters as the card writer and retain signed complex couplings
class AmplitudeSteering:
    # Resolve one trial without converting its differentiable coordinates to Python scalars
    def __init__(self, driver, path, parameters):
        self.values, physical, groups = driver._prepare_card_params(parameters, path=path if isinstance(path, dict) else str(path))
        self.names = physical["exact"]
        self.dependencies = {key: {physical["exact"].get(key, key)} for key in self.values}
        for group in (*groups[0], *groups[1], *groups[2]):
            sources = {physical["exact"].get(key, key) for key in group.keys}
            self.dependencies.update((key, sources) for key in group.targets)
        self.residues = {group.targets[0]: (group.amplitude, None) for group in groups[2]}
        for entry in groups[3].values():
            keys = [*entry["ordered_angle_keys"], *entry.get("cartesian_keys", [])]
            keys.extend(entry[key] for key in ("norm_key", "coupled_phase_key") if key in entry)
            sources = self.sources(keys)
            self.dependencies.update((key, sources) for key in (*entry["mag_keys"], *entry["phase_keys"]))
            self.residues.update((key, (value, entry.get("coupled_phase_key")))
                for key, value in zip(entry["mag_keys"], entry["amplitudes"], strict=True))

    # Identify the original optimizer coordinates represented by selected physical fields
    def sources(self, fields):
        return set().union(*(self.dependencies.get(key, {self.names.get(key, key)}) for key in fields))

    # Compute a production coupling including its common resonance phase exactly once
    def coupling(self, magnitude, angle, default, *, common=None, common_default=0.0):
        if magnitude in self.residues:
            value, included = self.residues[magnitude]
        else:
            mag = self.values.get(magnitude, default[0])
            # Width derived strengths are already contained in the native amplitude
            value = (1.0 if mag is None else mag) * phase(self.values.get(angle, default[1]))
            included = None
        if common is not None and included != common:
            value = value * phase(self.values.get(common, common_default))
        return value


# Construct unit production cards and their differentiable resonance coefficients
class ResonanceCoefficients:
    # Read immutable source cards once and derive the bank columns from normal active parameters
    def __init__(self, *, driver, path, model, resonances, pid, initial, decay=None, bounds=None, nodes=None):
        self.driver, self.path, self.model = driver, Path(path), model
        self.cards = {str(filename.relative_to(self.path)): load_json_file(filename, loader=pyjson5.load)
                      for filename in self.path.rglob("*.json")}
        self.decay = decay
        self.general = load_json_file(self.path / "GENERAL.json", loader=pyjson5.load)["PARAM_REGGE"]
        decays = load_json_file(self.path / "DECAYS.json", loader=pyjson5.load)
        self.resonances, self.axes, self.columns = {}, {}, []
        decoded = AmplitudeSteering(driver, path, initial)
        for name in resonances:
            card = load_json_file(self.path / "RES" / f"{name}.json", loader=pyjson5.load)
            res = card["PARAM_RES"]
            channels = decays[str(res["PDG"])]
            channel = next(key for key in channels if key.startswith("[") and sorted(pyjson5.loads(key)) == sorted(pid))
            production = driver._resonance_block(card, model)
            if model == "MP" and production["polarization"]["mode"] == "rho":
                raise ValueError("ampfit requires MP polarization.mode a_Jz or none, rho requires a sampling optimizer")
            pair = next(key for key, value in res["MODELS"][model].items() if value is production)
            hadronic = 22 not in pyjson5.loads(pair)
            form = res["MODELS"][model]["FF_prod"]
            if hadronic and form["type"] != "none" and (form["type"], form["norm"]) != ("gaussian", "pole"):
                raise ValueError(f"Unsupported differentiable resonance form factor for {name}")
            pole = any(f"RES|{name}:{model}:{field}" in decoded.values for field in ("mass", "width"))
            if pole and decay is None:
                raise ValueError("Pole reweighting requires decay normalization controls from ampfit.json")
            if pole and (not hadronic or res["MODELS"][model]["BW"] not in {"fixed-width", "kinematic-width"}):
                raise ValueError(f"Varying the pole of {name} requires hadronic production and a supported line shape")
            if pole and self.general["DECAY_BARRIERS"][model]:
                raise ValueError(f"Varying the pole of {name} requires decay barriers to be disabled")
            decay_form = channels[channel]["FF_decay"][model]
            if pole and decay_form["type"] != "none":
                raise ValueError(f"Varying the pole of {name} requires a differentiable decay form factor")
            self.resonances[name] = {"card": card, "decay": channels[channel], "channel": channel,
                                     "form": form, "hadronic": hadronic, "pole": pole,
                                     "analytic_prod": form["type"] == "gaussian" and form["norm"] == "pole"}
            entry = self.resonances[name]
            prefix = f"RES|{name}:{model}:"
            fields = (("FF_transfer.",) if entry["hadronic"] or entry["analytic_prod"]
                      else ("FF_transfer.", "FF_prod."))
            axes = {key: ChebyshevGrid(bound["lower"], bound["upper"], nodes)
                    for key, bound in sorted((bounds or {}).items())
                    if key.startswith(prefix) and key[len(prefix):].startswith(fields)}
            self.axes[name] = axes
            points = list(product(*(axis.nodes.tolist() for axis in axes.values())))
            self.columns.extend(dict(column, grid=index, parameters={
                **column["parameters"], **dict(zip(axes, point, strict=True))})
                for column in self._production(name, decoded) for index, point in enumerate(points))
        targets = {key for axes in self.axes.values() for key in axes}
        for name, entry in self.resonances.items():
            prefix = f"RES|{name}:"
            targets.add(f"{prefix}{model}:phi")
            if entry["pole"]:
                targets.update(f"{prefix}{model}:{field}" for field in ("mass", "width"))
            if entry["hadronic"]:
                targets.add(f"REGGE|omega.{model}")
            if entry["analytic_prod"]:
                targets.add(f"{prefix}{model}:FF_prod.Lambda2")
            decay_prefix = f"DECAY|{entry['card']['PARAM_RES']['PDG']}:{entry['channel']}:"
            targets.update(decay_prefix + field for field in ("BR", f"zeta.{model}"))
        for column in self.columns:
            if column["magnitude"] is not None:
                targets.update((column["magnitude"], column["angle"]))
            if "spin" in column:
                prefix = f"RES|{column['resonance']}:{model}:"
                targets.update(key for key in decoded.values if key.startswith(prefix)
                               and (driver._is_ajz_angle_key(key) or driver._is_ajz_phase_key(key)))
        self.used = decoded.sources(targets) & initial.keys()

    # Derive each physical production row and its unit coupling steering
    def _production(self, name, decoded):
        entry = self.resonances[name]
        res = entry["card"]["PARAM_RES"]
        block = self.driver._resonance_block(entry["card"], self.model)
        prefix = f"RES|{name}:{self.model}:"
        angular = any(key.startswith(prefix) and (self.driver._is_ajz_angle_key(key) or self.driver._is_ajz_phase_key(key))
                      for key in decoded.values)
        active = angular or any(key.startswith(prefix) and key[len(prefix):].startswith(("g", "helicity"))
                                for key in decoded.values)
        # Central form factors are applied after coherent screening and need no grid interpolation
        common = {prefix + "phi": 0.0}
        if entry["hadronic"] or entry["analytic_prod"]:
            common[prefix + "FF_prod"] = {"type": "none"}
        if not active:
            return [dict(resonance=name, parameters=common, magnitude=None)]
        if block["basis"].startswith("auto_"):
            common[prefix + "g[1]"] = 0.0
            if block["g"][0] is not None:
                common[prefix + "g[0]"] = 1.0
            coefficient = dict(magnitude=prefix + "g[0]", angle=prefix + "g[1]", default=block["g"])
            columns = [dict(resonance=name, parameters=common, **coefficient)]
        else:
            field = block["basis"]
            if field not in {"g_ls", "helicity"}:
                raise ValueError(f"Unsupported amplitude production basis for {name}: {block['basis']}")
            rows = block[field]
            zero = {prefix + angular_row_name(field, row) + suffix: 0.0 for row in rows for suffix in ("@MAG", "@PHASE")}
            columns = [dict(
                resonance=name, parameters={**common, **zero, prefix + angular_row_name(field, row) + "@MAG": 1.0},
                magnitude=f"{prefix}{field}[{index},2]", angle=f"{prefix}{field}[{index},3]", default=row[2:4]
            ) for index, row in enumerate(rows)]
        if not angular:
            return columns
        spin = res["spinX2"] // 2
        sectors = mp_spin_sectors(spin, decoded.values.get("REGGE|MP_FRAME", self.general["MP_FRAME"]))
        # Tensor products of coupling rows and pure projectors span the steered amplitude
        states = [(m, m, 0) for m in sectors]
        states += [(m, n, q) for m, n in combinations(sectors, 2) for q in (0, 1)]
        return [dict(column, spin=(m, n, q), parameters={**column["parameters"], prefix + "polarization.mode": "a_Jz",
            prefix + "polarization.a_Jz": [
                [-k, float(k in (m, n)) / math.sqrt((1 if k == 0 else 2) * (1 if m == n else 2)),
                 q * math.pi / 2 if k == n else 0.0] for k in range(spin, -1, -1)]
        }) for column, (m, n, q) in product(columns, states)]

    # Compute the scalar factors which multiply all production rows of one resonance
    def _factors(self, name, values, mass2, daughters):
        entry = self.resonances[name]
        res, decay = entry["card"]["PARAM_RES"], entry["decay"]
        pole = res["MODELS"][self.model]
        model_prefix = f"RES|{name}:{self.model}:"
        mass, width = values.get(model_prefix + "mass", pole["mass"]), values.get(model_prefix + "width", pole["width"])
        form = entry["form"]
        factor = torch.ones_like(mass2, dtype=torch.complex128)
        if entry["analytic_prod"]:
            scale = values.get(model_prefix + "FF_prod.Lambda2", form["Lambda2"])
            factor = factor * torch.exp(-((mass2 - mass**2) / scale)**2)
        if entry["pole"]:
            imaginary = torch.sqrt(mass2) if pole["BW"] == "kinematic-width" else mass
            imaginary0 = torch.sqrt(mass2) if pole["BW"] == "kinematic-width" else pole["mass"]
            factor = factor * (mass2 - pole["mass"]**2 + 1j * imaginary0 * pole["width"]) / (mass2 - mass**2 + 1j * imaginary * width)
            # The fixed spin multiplicity and two body symmetry factor cancel in the decay normalization ratio
            phase_space = decay_phase_space(mass, width, daughters, self.decay)
            phase_space0 = float(decay_phase_space(pole["mass"], pole["width"], daughters, self.decay))
            factor = factor * array.sqrt(width / pole["width"] * phase_space0 / phase_space)
        decay_prefix = f"DECAY|{res['PDG']}:{entry['channel']}:"
        factor = factor * array.sqrt(values.get(decay_prefix + "BR", decay["BR"]) / decay["BR"])
        phase0 = decay["zeta"][self.model]
        factor = factor * phase(values.get(decay_prefix + f"zeta.{self.model}", phase0) - phase0)
        if entry["hadronic"]:
            omega0 = self.general["omega"][self.model]
            factor = factor * (self.general["s0"] / mass2)**(values.get(f"REGGE|omega.{self.model}", omega0) - omega0)
        return factor

    # Interpolate each resonance after removing its known external transfer dependence
    def _grid(self, name, values, mass2, transfer):
        weights = torch.ones((len(mass2), 1), dtype=mass2.dtype, device=mass2.device)
        for key, axis in self.axes[name].items():
            cardinal = axis(values[key]).to(mass2.device).expand(len(mass2), -1)
            if key.endswith(("FF_transfer.Lambda2", "FF_transfer.LambdaInv2")):
                form = self.resonances[name]["card"]["PARAM_RES"]["MODELS"][self.model]["FF_transfer"]
                if (form["type"], form["norm"]) == ("power", "zero"):
                    cardinal = transfer_weights(axis, values[key], transfer, form["n"], key.endswith("LambdaInv2"))
            weights = (weights[:, :, None] * cardinal[:, None, :]).flatten(1)
        return weights

    # Evaluate all complex production coefficients with shared phases and active parameter transforms
    def evaluate(self, parameters, mass2, daughters, transfer, basis=None):
        steering = parameters if isinstance(parameters, AmplitudeSteering) else AmplitudeSteering(self.driver, self.path, parameters)
        factors = {name: self._factors(name, steering.values, mass2, daughters) for name in self.resonances}
        grids = {name: self._grid(name, steering.values, mass2, transfer) for name in self.axes}
        spins, coefficients = {}, []
        for column in self.columns:
            if basis is not None and column["grid"] != 0:
                continue
            name = column["resonance"]
            res = self.resonances[name]["card"]
            prefix = f"RES|{name}:{self.model}:"
            common = prefix + "phi"
            default_phase = res["PARAM_RES"]["MODELS"][self.model]["phi"]
            coefficient = (phase(steering.values.get(common, default_phase)) if column["magnitude"] is None else
                           steering.coupling(column["magnitude"], column["angle"], column["default"],
                                             common=common, common_default=default_phase))
            if "spin" in column:
                if name not in spins:
                    geometry = "projective" if any(key.startswith(prefix) and key.endswith("@AJZP") for key in steering.values) else "sphere"
                    spins[name] = self.driver._ajz_sectors(res=name, j=res["PARAM_RES"]["spinX2"] // 2,
                                                        geometry=geometry, param={"REGGE|MP_FRAME": self.general["MP_FRAME"], **steering.values}, json_data=res)
                    spins[name] = {m: array.asarray(value, dtype=complex) for m, value in spins[name].items()}
                state = spins[name]
                m, n, quadrature = column["spin"]
                if m == n:
                    weight = (state[m] * state[m].conj()).real
                    for other in state:
                        if other == m:
                            continue
                        pair = state[min(m, other)] * state[max(m, other)].conj()
                        weight = weight - pair.real + pair.imag
                else:
                    pair = state[m] * state[n].conj()
                    weight = 2 * (-pair.imag if quadrature else pair.real)
                coefficient = coefficient * weight
            coefficients.append(torch.as_tensor(coefficient, device=mass2.device))
        if basis is None:
            return torch.stack([coefficient * factors[column["resonance"]] * grids[column["resonance"]][:, column["grid"]]
                                for coefficient, column in zip(coefficients, self.columns, strict=True)], dim=1)
        # Contract production rows before applying their shared event dependent factors
        amplitudes, offset, row = [], 0, 0
        for name in self.resonances:
            count = sum(column["resonance"] == name for column in self.columns)
            grid = grids[name]
            block = basis[offset:offset + count].reshape(-1, grid.shape[1], *basis.shape[1:])
            rows = count // grid.shape[1]
            couplings = torch.stack(coefficients[row:row + rows]).to(basis)
            projected = torch.einsum("rknh,r->knh", block, couplings)
            amplitudes.append(factors[name][:, None] * torch.einsum("knh,nk->nh", projected, grid.to(basis)))
            offset += count
            row += rows
        return torch.stack(amplitudes).sum(0)


# Compute the logarithm of one pole normalized meson form factor
def offshell_log(ff_type, parameters, virtuality):
    if ff_type == "exp":
        return -parameters["b"] * virtuality
    if ff_type == "power":
        return -parameters["n"] * torch.log1p(virtuality / parameters["Lambda2"])
    if ff_type == "logexp":
        y = torch.log1p(virtuality / parameters["Lambda2"])
        return -parameters["b"] * parameters["Lambda2"] * y * (1 + 0.5 * y)
    a = parameters["a"]
    return -parameters["b"] * (torch.sqrt(virtuality + a * a) - a)


# Derive continuum coupling products and screened form factor interpolation from source cards
class ContinuumCoefficients:
    # Resolve shared JSON parameters and enumerate the same exchange pairs as the generator
    def __init__(self, *, driver, path, model, pid, initial, bounds, nodes):
        self.driver, self.path = driver, Path(path)
        filename = self.path / f"CON_{model}.json"
        reader = jsonref.JsonReader(pyjson5.load)
        source = reader.read(filename)
        general = load_json_file(self.path / "GENERAL.json", loader=pyjson5.load)["PARAM_REGGE"]
        self.channel, _, self.pairs = _continuum_entry_for_pair(general["PARAM_CON"][model], tuple(pid))
        if len(pid) != 2 or pid[0] != -pid[1] or any(22 in pair for pair in self.pairs):
            raise ValueError("Differentiable continuum banks require a charged meson pair and hadronic exchanges")
        decoded = AmplitudeSteering(driver, path, initial)
        origins = {}
        for key in decoded.values:
            if not key.startswith(f"CON_{model}|"):
                continue
            _, exchange, channel, field = key.replace("|", ":", 1).split(":")
            field = field.replace("LambdaInv2", "Lambda2")
            target, parts = reader.origin(filename, [exchange, *channel.split("/"), *jsonref.field_parts(field)])
            origins[target, tuple(parts)] = key

        # Match each vertex field to the active coordinate which edits its shared source
        def parameter(exchange, family, field):
            target, parts = reader.origin(filename, [exchange, family, *jsonref.field_parts(field)])
            return origins.get((target, tuple(parts)))

        exchanges = list(dict.fromkeys(str(exchange) for pair in self.pairs for exchange in pair))
        self.vertices, self.forms, self.transfers, fields, shared, disable_veto = {}, {}, {}, {}, {}, {}
        unit = {}
        veto_keys = set()
        for exchange in exchanges:
            family, block = _production_entry_for_pair(source[exchange], tuple(pid), card_name=filename.name)
            prefix = f"CON_{model}|{exchange}:{family}"
            vertex = driver.production_sector_block(source[exchange], family + "/opposite")
            base = prefix + "/opposite:"
            if "g" in vertex:
                rows = [(base + "g[0]", base + "g[1]", vertex["g"])]
            else:
                field = {"crossed_helicity": "helicity", "crossed_ls": "g_ls"}.get(vertex.get("basis"))
                if field is None or not vertex.get(field):
                    raise ValueError("Differentiable meson continuum requires an explicit coupling table")
                position, _ = driver._angular_columns(base + field)
                rows = [(f"{base}{field}[{index},{position}]", f"{base}{field}[{index},{position + 1}]", row[-2:])
                        for index, row in enumerate(vertex[field])]
            if any(default[0] is None or default[0] < 0 for _, _, default in rows):
                raise ValueError("Continuum source couplings must be nonnegative")
            self.vertices[exchange] = [row for row in rows if row[2][0] > 0 or row[0] in decoded.values]
            if not self.vertices[exchange]:
                raise ValueError("The continuum vertex must contain an active coupling")
            unit[exchange] = {key: 0.0 for magnitude, angle, _ in rows for key in (magnitude, angle)}
            offshell, transfer = block["FF_offshell"], block["FF_transfer"]
            if offshell["type"] not in {"exp", "power", "orear", "logexp"} or offshell["norm"] != "pole" or (transfer["type"], transfer["norm"]) != ("power", "zero"):
                raise ValueError("Continuum grids require exp, power, orear or logexp pole off shell and power transfer form factors")
            ff_fields = {field: parameter(exchange, family, f"FF_offshell.{field}")
                         for field in offshell if field not in {"type", "norm"}}
            self.forms[exchange] = (offshell, ff_fields)
            self.transfers[exchange] = (transfer, parameter(exchange, family, "FF_transfer.Lambda2"))
            for field in (*("FF_offshell." + name for name in ff_fields), "FF_transfer.Lambda2"):
                key = parameter(exchange, family, field)
                if key in bounds:
                    fields[key] = bounds[key]
                    shared.setdefault(field, set()).add(key)
            for field in ("pveto.M0", "pveto.c"):
                veto_keys.add((field, parameter(exchange, family, field)))
            disable_veto[prefix + ":pveto.active"] = False
        if len(fields) > 3 or any(len(keys) > 1 for keys in shared.values()) or len(veto_keys) > 2:
            raise ValueError("Screened continuum interpolation requires shared off shell, transfer and veto parameters")
        if not isinstance(nodes, int) or isinstance(nodes, bool) or nodes < 2:
            raise ValueError("Continuum interpolation needs at least two nodes per active form factor")
        self.veto = block["pveto"]
        self.veto_keys = dict(veto_keys)
        self.axes = {key: ChebyshevGrid(bound["lower"], bound["upper"], nodes) for key, bound in sorted(fields.items())}
        points = list(product(*(axis.nodes.tolist() for axis in self.axes.values())))
        self.grid = {key: torch.tensor([point[index] for point in points], dtype=torch.float64)
                     for index, key in enumerate(self.axes)}
        # Equal exchanges share one table, so mixed rows use the quadratic polarization identity
        self.terms, self.columns = [], []
        for pair in self.pairs:
            first, second = map(str, pair)
            indices = range(len(self.vertices[first]))
            rows = (combinations_with_replacement(indices, 2) if first == second
                    else product(indices, range(len(self.vertices[second]))))
            for i, j in rows:
                self.terms.append((first, second, i, j))
                parameters = {**disable_veto, **unit[first], **unit[second],
                              self.vertices[first][i][0]: 1.0, self.vertices[second][j][0]: 1.0}
                self.columns.extend(dict(exchange_pair=pair, parameters={**parameters, **dict(zip(self.axes, point, strict=True))})
                                    for point in points)
        targets = {*self.axes, *(key for rows in self.vertices.values() for magnitude, angle, _ in rows
                                for key in (magnitude, angle))}
        if self.veto["active"]:
            targets.update(key for key in self.veto_keys.values() if key is not None)
        self.used = decoded.sources(targets) & initial.keys()

    # Contract physical exchange couplings with the fixed interpolation polynomials and central veto
    def evaluate(self, parameters, mass2, virtuality, transfer, basis=None):
        steering = parameters if isinstance(parameters, AmplitudeSteering) else AmplitudeSteering(self.driver, self.path, parameters)
        weights = torch.ones((len(mass2), 1), dtype=mass2.dtype, device=mass2.device)
        for key, axis in self.axes.items():
            cardinal = axis(steering.values[key]).to(mass2.device).expand(len(mass2), -1)
            weights = (weights[:, :, None] * cardinal[:, None, :]).flatten(1)
        corrections = {}
        for exchange, (form, keys) in self.forms.items():
            current = {field: steering.values.get(key, form[field]) for field, key in keys.items()}
            nodes = {field: self.grid[key].to(mass2.device) if key in self.grid else current[field]
                     for field, key in keys.items()}
            x = virtuality[:, None] / 2
            correction = torch.exp(offshell_log(form["type"], current, x) - offshell_log(form["type"], nodes, x))
            # Extract the external transfer factor on each beam leg before interpolating screening
            form, key = self.transfers[exchange]
            correction = correction[:, :, None].expand(-1, -1, 2)
            if key in self.grid:
                value, grid = steering.values[key], self.grid[key].to(mass2.device)
                if not key.endswith("LambdaInv2"):
                    value, grid = 1 / value, grid.reciprocal()
                ratio = (1 - grid[None, :, None] * transfer[:, None, :]) / (1 - value * transfer[:, None, :])
                correction = correction * ratio.pow(form["n"])
            corrections[exchange] = correction
        vertices = {exchange: [steering.coupling(magnitude, angle, default) for magnitude, angle, default in rows]
                    for exchange, rows in self.vertices.items()}
        couplings = []
        for first, second, i, j in self.terms:
            value = vertices[first][i] * vertices[second][j]
            if first == second and i == j:
                value = value - vertices[first][i] * sum(c for k, c in enumerate(vertices[first]) if k != i)
            couplings.append(value)
        veto = torch.ones_like(mass2)
        if self.veto["active"]:
            mass = torch.sqrt(mass2)
            mass0 = steering.values.get(self.veto_keys["pveto.M0"], self.veto["M0"])
            power = steering.values.get(self.veto_keys["pveto.c"], self.veto["c"])
            veto = torch.where(mass <= mass0, torch.ones_like(mass), (mass0 / mass)**power)
        if basis is None:
            return torch.cat([veto[:, None] * weights * corrections[first][:, :, 0] * corrections[second][:, :, 1] * value
                              for (first, second, _, _), value in zip(self.terms, couplings, strict=True)], dim=1)
        blocks = basis.reshape(len(self.terms), -1, *basis.shape[1:])
        return veto[:, None] * sum(value * torch.einsum(
            "knh,nk->nh", block, (weights * corrections[first][:, :, 0] * corrections[second][:, :, 1]).to(basis))
            for block, (first, second, _, _), value in zip(blocks, self.terms, couplings, strict=True))
