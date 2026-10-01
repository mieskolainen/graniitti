#!/usr/bin/env python3
# Derive continuum couplings from Donnachie-Landshoff fits
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

# sigma(h+/- p) = P*s^eps + (C+ -/+ C-)*s^-eta
# SOFT fixes g_p, DL fixes g_h/g_p = C_hp/C_pp
# Pole convention: ||H||^2/(2*s_h+1) = g_h^2, GP normalized at m=0
# Preserve relative phases and rescale each allowed charge sector
# [REFERENCE: arXiv:hep-ph/9209205, arXiv:1804.04706, Eq. (3.29)]

from __future__ import annotations

import argparse
import cmath
import copy
import json
import math
import os
import re
import stat
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any
from urllib.parse import unquote

# Resolve the tools package for direct execution
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from develop.tools.lib import common, soft_exchange
from develop.tools.lib import pole as spinmath


# Card parsing and scalar updates
@dataclass(frozen=True)
class ScalarUpdate:
    path: tuple[Any, ...]
    label: str
    value: Any


class CardDocument:
    """JSON5 values and scalar spans for exact card edits"""
    # Parse card values and scalar spans
    def __init__(self, source: str, label: str = "<card>"):
        self.source, self.label, self.pos = source, label, 0
        self.spans: dict[tuple[Any, ...], tuple[int, int]] = {}
        self.data = self._value(())
        self._space()
        if self.pos != len(source):
            self._error("unexpected trailing content")

    # Report a parser error with its source line
    def _error(self, text: str):
        line = self.source.count("\n", 0, self.pos) + 1
        raise ValueError(f"{self.label}:{line}: {text}")

    # Skip whitespace and comments
    def _space(self):
        s = self.source
        while self.pos < len(s):
            if s[self.pos].isspace() or s[self.pos] == '\ufeff':
                self.pos += 1
            elif s.startswith('//', self.pos):
                end = s.find('\n', self.pos + 2)
                self.pos = len(s) if end < 0 else end + 1
            elif s.startswith('/*', self.pos):
                end = s.find('*/', self.pos + 2)
                if end < 0:
                    self._error("unterminated comment")
                self.pos = end + 2
            else:
                break

    # Decode a JSON5 string
    def _string(self) -> str:
        s, quote = self.source, self.source[self.pos]
        self.pos += 1
        out = []
        escapes = {'b': '\b', 'f': '\f', 'n': '\n', 'r': '\r', 't': '\t', 'v': '\v', '0': '\0'}
        while self.pos < len(s):
            char = s[self.pos]
            self.pos += 1
            if char == quote:
                # Join valid UTF-16 surrogate pairs from JSON unicode escapes
                text = ''.join(out)
                return text.encode('utf-16', 'surrogatepass').decode('utf-16', 'surrogatepass')
            if char in '\r\n':
                self._error("unescaped newline in string")
            if char != '\\':
                out.append(char)
                continue
            if self.pos >= len(s):
                self._error("unterminated string escape")
            char = s[self.pos]
            self.pos += 1
            if char in ('u', 'x'):
                size = 4 if char == 'u' else 2
                digits = s[self.pos:self.pos+size]
                if len(digits) != size or not re.fullmatch(r'[0-9a-fA-F]+', digits):
                    self._error("invalid hexadecimal string escape")
                out.append(chr(int(digits, 16)))
                self.pos += size
            elif char == '\r':
                if self.pos < len(s) and s[self.pos] == '\n':
                    self.pos += 1
            elif char == '\n':
                pass
            elif char == '0' and self.pos < len(s) and s[self.pos].isdigit():
                self._error("octal string escapes are not supported")
            else:
                out.append(escapes.get(char, char))
        self._error("unterminated string")

    # Parse a value and retain scalar spans
    def _value(self, path: tuple[Any, ...]):
        self._space()
        if self.pos >= len(self.source):
            self._error("expected a value")
        start, char = self.pos, self.source[self.pos]
        if char in '{[':
            self.pos += 1
            obj = {} if char == '{' else []
            close = '}' if char == '{' else ']'
            self._space()
            while self.pos < len(self.source) and self.source[self.pos] != close:
                if char == '{':
                    if self.source[self.pos] in "\"'":
                        key = self._string()
                    else:
                        match = re.match(r'[A-Za-z_$][\w$]*', self.source[self.pos:])
                        if not match:
                            self._error("expected a quoted key or identifier")
                        key = match.group()
                        self.pos += len(key)
                    if key in obj:
                        self._error(f"duplicate object key {key!r}")
                    self._space()
                    if self.pos >= len(self.source) or self.source[self.pos] != ':':
                        self._error("expected ':'")
                    self.pos += 1
                    obj[key] = self._value((*path, key))
                else:
                    obj.append(self._value((*path, len(obj))))
                self._space()
                if self.pos >= len(self.source):
                    self._error(f"expected '{close}'")
                if self.source[self.pos] == close:
                    break
                if self.source[self.pos] != ',':
                    self._error("expected ','")
                self.pos += 1
                self._space()
            if self.pos >= len(self.source):
                self._error(f"expected '{close}'")
            self.pos += 1
            return obj
        if char in "\"'":
            value = self._string()
        else:
            match = re.match(r'(?:true|false|null|[+-]?(?:Infinity|NaN|0[xX][0-9a-fA-F]+|(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?))',
                             self.source[self.pos:])
            if not match:
                self._error("invalid value")
            token = match.group()
            self.pos += len(token)
            if token in ('true', 'false', 'null'):
                value = {'true': True, 'false': False, 'null': None}[token]
            elif 'x' in token.lower():
                value = int(token, 16)
            elif re.fullmatch(r'[+-]?\d+', token):
                value = int(token)
            else:
                value = float(token.replace('Infinity', 'inf').replace('NaN', 'nan'))
        self.spans[path] = (start, self.pos)
        return value


# Read a value at a card path
def _at(value: Any, path: tuple[Any, ...]):
    for key in path:
        value = value[int(key)] if isinstance(value, list) else value[key]
    return value


# Assign a value at a card path
def _set(value: Any, path: tuple[Any, ...], replacement: Any):
    if not path:
        raise ValueError("cannot replace an entire card")
    parent = _at(value, path[:-1])
    parent[int(path[-1]) if isinstance(parent, list) else path[-1]] = replacement


class CardReader:
    """Resolve local references within selected cards"""
    # Set allowed card paths and an empty cache
    def __init__(self, allowed):
        self.allowed = {Path(path).resolve() for path in allowed}
        self.docs: dict[Path, CardDocument] = {}

    # Cache parsed cards within the allowed paths
    def document(self, path):
        path = Path(path).resolve()
        if path not in self.allowed:
            raise ValueError(f"External input is not permitted: {path}. Inline it in GENERAL.json or a selected CON card.")
        if path not in self.docs:
            with path.open('r', encoding='utf-8', newline='') as stream:
                self.docs[path] = CardDocument(stream.read(), str(path))
        return self.docs[path]

    # Resolve a local card reference
    def _reference(self, card: Path, value: Any):
        ref = value['$ref']
        if not isinstance(ref, str) or set(value) != {'$ref'}:
            raise ValueError(f"{card}: $ref must be a standalone string reference")
        name, _, pointer = ref.partition('#')
        if '://' in name or name.startswith('//'):
            raise ValueError("Network references are not permitted")
        target = (card.parent / name).resolve() if name else card
        self.document(target)  # Enforce the complete input allowlist
        pointer = unquote(pointer)
        if pointer and not pointer.startswith('/'):
            raise ValueError(f"{card}: $ref fragment must be a JSON pointer")
        parts = tuple(x.replace('~1', '/').replace('~0', '~') for x in pointer.split('/')[1:]) if pointer else ()
        return target, parts

    # Locate the source of a referenced value
    def origin(self, card, path, seen=frozenset()):
        card = Path(card).resolve()
        value = self.document(card).data
        walked = ()
        remaining = tuple(path)
        while True:
            if isinstance(value, dict) and '$ref' in value:
                key = (card, walked)
                if key in seen:
                    raise ValueError(f"Cyclic $ref in {card}")
                target, parts = self._reference(card, value)
                return self.origin(target, parts + remaining, seen | {key})
            if not remaining:
                return card, walked
            key, *tail = remaining
            key = int(key) if isinstance(value, list) else key
            value = value[key]
            walked += (key,)
            remaining = tuple(tail)

    # Read a card with references expanded
    def read(self, card):
        # Expand references recursively
        def expand(file, parts, stack):
            file, parts = self.origin(file, parts)
            key = (file, parts)
            if key in stack:
                raise ValueError(f"Cyclic $ref in {file}")
            value = _at(self.document(file).data, parts)
            stack = stack | {key}
            if isinstance(value, dict):
                return {name: expand(file, (*parts, name), stack) for name in value}
            if isinstance(value, list):
                return [expand(file, (*parts, i), stack) for i in range(len(value))]
            return value
        return expand(Path(card).resolve(), (), frozenset())


# Read a card through the current derivation reader
def read_card(path, reader=None):
    path = Path(path).resolve()
    return (reader or CardReader({path})).read(path)


# Resolve writable coupling locations while retaining read-only input dependencies
def _prepare_updates(updates_by_card, reader=None):
    reader = reader or CardReader(updates_by_card)
    writable = {Path(path).resolve() for path in updates_by_card}
    desired, labels, assigned = {}, {}, {}
    for card, updates in updates_by_card.items():
        for update in updates:
            target, parts = reader.origin(card, update.path)
            if target not in writable:
                raise ValueError(f"A coupling points outside the selected target CON cards: {target}. GENERAL.json is read-only.")
            key = (target, parts)
            if key in assigned:
                if assigned[key] != update.value:
                    raise ValueError(f"Conflicting updates for shared parameter {target}:{parts}")
                continue
            if isinstance(update.value, float) and not math.isfinite(update.value):
                raise ValueError(f"Nonfinite output coupling {update.label}")
            old = _at(reader.document(target).data, parts)
            if isinstance(old, (dict, list)):
                raise ValueError(f"Update is not a scalar: {update.label}")
            assigned[key], labels[key] = update.value, update.label
            desired.setdefault(target, {})[parts] = update.value
    outputs, rows, sources = {}, [], {}
    for card, changes in desired.items():
        document = reader.document(card)
        sources[card] = document.source
        edits = []
        for path, value in changes.items():
            old = _at(document.data, path)
            if old == value:
                continue
            start, end = document.spans[path]
            edits.append((start, end, json.dumps(value, allow_nan=False, ensure_ascii=False)))
            rows.append((labels[card, path], old, value))
        if edits:
            source = document.source
            for start, end, replacement in sorted(edits, reverse=True):
                source = source[:start] + replacement + source[end:]
            CardDocument(source, str(card))
            outputs[card] = source
    return outputs, rows, sources


# Read card text with original line endings
def _read_exact(path):
    with Path(path).open('r', encoding='utf-8', newline='') as stream:
        return stream.read()


# Apply confirmed scalar updates through the same reader used for their derivation
def push_json5_updates(updates_by_card, *, confirm=None, dry_run=False, reader=None):
    reader = reader or CardReader(updates_by_card)
    outputs, rows, sources = _prepare_updates(updates_by_card, reader)
    if not rows:
        print('No parameter changes to push')
        return False
    print(common.table(['parameter', 'old', 'new'],
                       [(name, f'{old:.12g}' if isinstance(old, (float, int)) else old,
                         f'{new:.12g}' if isinstance(new, (float, int)) else new) for name, old, new in rows],
                       right_align={1, 2}))
    if dry_run:
        print('Dry run: no files written')
        return False
    if confirm is None:
        try:
            accepted = input(f'Update {len(outputs)} target card(s)? [y/N] ').strip().lower() in {'y', 'yes'}
        except (EOFError, KeyboardInterrupt):
            accepted = False
    else:
        accepted = confirm(rows)
    if not accepted:
        print('Push cancelled')
        return False
    # Recheck every read dependency, including the read-only GENERAL card
    sources.update({p: doc.source for p, doc in reader.docs.items()})
    for card, source in sources.items():
        if _read_exact(card) != source:
            raise ValueError(f'Card changed after preview: {card}')
    staged, replaced, backups = {}, [], []
    try:
        for card, text in outputs.items():
            backup = card.with_name(card.name + '.bak')
            suffix = 0
            while backup.exists():
                suffix += 1
                backup = card.with_name(card.name + f'.bak.{suffix}')
            with backup.open('x', encoding='utf-8', newline='') as stream:
                stream.write(sources[card])
            backups.append(backup.resolve())
            descriptor, name = tempfile.mkstemp(prefix='.' + card.name + '.', suffix='.tmp', dir=card.parent)
            with os.fdopen(descriptor, 'w', encoding='utf-8', newline='') as stream:
                stream.write(text)
                stream.flush()
                os.fsync(stream.fileno())
            os.chmod(name, stat.S_IMODE(card.stat().st_mode))
            staged[card] = Path(name)
        for card, temp in staged.items():
            if _read_exact(card) != sources[card]:
                raise ValueError(f'Card changed during write preparation: {card}')
            os.replace(temp, card)
            replaced.append(card)
    except BaseException:
        # Per-file replacement is atomic. Restore earlier replacements on a
        # detected failure; a process/OS crash is not a multi-file transaction
        for card in reversed(replaced):
            descriptor, name = tempfile.mkstemp(prefix='.' + card.name + '.', dir=card.parent)
            with os.fdopen(descriptor, 'w', encoding='utf-8', newline='') as stream:
                stream.write(sources[card])
            os.chmod(name, stat.S_IMODE(card.stat().st_mode))
            os.replace(name, card)
        raise
    finally:
        for temp in staged.values():
            temp.unlink(missing_ok=True)
    print(f'Push completed: updated {len(rows)} values in {len(outputs)} card files')
    print('Backups:')
    for backup in backups:
        print(f'  {backup}')
    return True


# Standard finite-spin aliases from the continuum models. GP uses pole_spin
# in GENERAL.json. Natural parity P=(-1)^J is explicit, not fitted here
_FIXED_EXCHANGE_SPINS = {991: 0, 993: 1, 995: 2, 9915: 2, 9925: 2,
                         9933: 1, 9943: 1, 9993: 1}
_ANALYTIC_EXCHANGE_PDGS = {990, 9910, 9920, 9930, 9940, 9990}


# Load built-ins and GENERAL.PARAM_PDG
def load_particles(tune_dir: Path, reader=None) -> dict[int, spinmath.Particle]:
    general = read_card(Path(tune_dir) / "GENERAL.json", reader)
    config = soft_exchange.load(general)
    particles = spinmath.particles()
    for row in config.regge_exchanges:
        crossing = int(config.definitions[row.soft_exchange]["crossing"])
        for pdg in row.pdg:
            if pdg in _FIXED_EXCHANGE_SPINS:
                spin = _FIXED_EXCHANGE_SPINS[pdg]
            elif pdg in _ANALYTIC_EXCHANGE_PDGS:
                spin = row.pole_spin
            else:
                continue  # A nonstandard alias must have explicit inline metadata
            particles[pdg] = spinmath.Particle(pdg, f"{row.soft_exchange}:{pdg}", 2 * spin,
                                      1 if spin % 2 == 0 else -1, crossing)
    # Optional customization stays inside the one permitted GENERAL input
    inline = general.get("PARAM_PDG", {})
    if not isinstance(inline, dict):
        raise ValueError("GENERAL.PARAM_PDG must be an object when supplied")
    for name, row in inline.items():
        if not isinstance(row, dict):
            raise ValueError(f"GENERAL.PARAM_PDG.{name} must be an object")
        pdg = soft_exchange.exact_integer(row.get("PDG"), f"PARAM_PDG.{name}.PDG")
        spin2 = soft_exchange.exact_integer(row.get("spinX2"), f"PARAM_PDG.{name}.spinX2")
        parity = soft_exchange.exact_integer(row.get("P"), f"PARAM_PDG.{name}.P")
        cparity = soft_exchange.exact_integer(row.get("C", 0), f"PARAM_PDG.{name}.C")
        if parity not in {-1,1} or cparity not in {-1,0,1}:
            raise ValueError(f"Invalid parity in GENERAL.PARAM_PDG.{name}")
        if spin2 < 0:
            matched = [r for r in config.regge_exchanges if pdg in r.pdg]
            if len(matched) != 1:
                raise ValueError(f"Analytic PDG {pdg} requires one GENERAL Regge mapping")
            spin2 = 2 * matched[0].pole_spin
        particles[pdg] = spinmath.Particle(pdg, str(row.get("name", name)), spin2, parity, cparity)
        if pdg > 0 and cparity == 0:
            particles[-pdg] = replace(particles[pdg], pdg=-pdg,
                                       parity=parity if spin2 % 2 == 0 else -parity)
    return particles


# Load Regge pole metadata
def load_regge(tune_dir: Path, reader=None) -> tuple[None, spinmath.ReggeTable]:
    config = soft_exchange.load(read_card(Path(tune_dir) / "GENERAL.json", reader))
    return None, spinmath.ReggeTable(
        groups=[list(row.pdg) for row in config.regge_exchanges],
        tau=[int(config.definitions[row.soft_exchange]["tau"]) for row in config.regge_exchanges],
        pole_spin=[row.pole_spin for row in config.regge_exchanges],
    )


# DL couplings and continuum cards
@dataclass(frozen=True)
class FinalState:
    """Continuum final state"""

    key: str
    label: str
    pdg: tuple[int | None, int | None]
    direct_hadron: str | None


@dataclass(frozen=True)
class Exchange:
    """DL exchange"""

    label: str
    key: str


@dataclass(frozen=True)
class ConfiguredExchange:
    """One mapped SOFT exchange and its physical proton projection"""

    exchange: Exchange
    soft_exchange: str
    trajectory_mode: str
    a0: float
    ap: float
    B: float
    beam_residue_t0: float
    beam_residue_source: str
    tau: int
    eta_mode: str


MB_TO_GEV2 = 2.56819
# [REFERENCE: DL, arXiv:hep-ph/9209205]
DL_EPSILON = 0.0808
DL_ETA = 0.4525
ROOT = next((p for p in Path(__file__).resolve().parents
             if (p / "modeldata" / "TUNE0" / "GENERAL.json").is_file()), Path.cwd())
GENERAL_CARD = ROOT / "modeldata" / "TUNE0" / "GENERAL.json"

POMERON = Exchange("Pomeron", "P")
F2_REGGEON = Exchange("C-even Reggeon", "C+")
RHO_REGGEON = Exchange("C-odd Reggeon", "C-")
EXCHANGES = (POMERON, F2_REGGEON, RHO_REGGEON)
ODDERON = Exchange("Odderon", "O")
CONFIGURED_EXCHANGES = (*EXCHANGES, ODDERON)

MODEL_EXCHANGE_PDGS = {
    "MP": {
        "P": (991, 993, 995),
        "C+": (9915, 9925),
        "C-": (9933, 9943),
    },
    "XP": {
        "P": (991, 993, 995),
        "C+": (9915, 9925),
        "C-": (9933, 9943),
    },
    "GP": {
        "P": (990,),
        "C+": (9910, 9920),
        "C-": (9930, 9940),
    },
}

ODD_SPIN_EXCHANGE_PDGS = {993, 9933, 9943, 9930, 9940}

TOOL_VERSION = "3.0.0"
COUPLING_PRECISION = 9
DEFAULT_VERTEX_NORMALIZATION = "pole-rms"
DEFAULT_SHAPE = "preserve"
VERTEX_NORMALIZATIONS = (DEFAULT_VERTEX_NORMALIZATION,)

FINAL_STATES = (
    FinalState("p", "p pbar", (2212, -2212), "p"),
    FinalState("n", "n nbar", (2112, -2112), "p"),
    FinalState("pi", "pi+ pi-", (211, -211), "pi"),
    FinalState("pi0", "pi0 pi0", (111, 111), "pi"),
    FinalState("K", "K+ K-", (321, -321), "K"),
    FinalState("K0", "K0 K0bar", (311, -311), "K0"),
    FinalState("rho", "rho0 rho0", (113, 113), None),
    FinalState("phi", "phi phi", (333, 333), None),
    FinalState("generic", "generic", (None, None), None),
)

# DL total-cross-section coefficients in mb:
# sigma = P * s^eps + Y * s^-eta
#
# For p, "plus" denotes pp and "minus" denotes pbar p
# For pi and K, "plus" denotes h+ p and "minus" denotes h- p
# For K0, use K0 p = K+ n and K0bar p = K- n by isospin
DL_MB = {
    # [REFERENCE: DL, arXiv:hep-ph/9209205]
    "p": {
        "pdg_pair": "2212,-2212",
        "plus": "p p",
        "minus": "pbar p",
        "P": 21.70,
        "Y_plus": 56.08,
        "Y_minus": 98.39,
    },
    # [REFERENCE: DL, arXiv:hep-ph/9209205]
    "pi": {
        "pdg_pair": "211,-211",
        "plus": "pi+ p",
        "minus": "pi- p",
        "P": 13.63,
        "Y_plus": 27.56,
        "Y_minus": 36.02,
    },
    # [REFERENCE: LNS, arXiv:1804.04706]
    "K": {
        "pdg_pair": "321,-321",
        "plus": "K+ p",
        "minus": "K- p",
        "P": 11.93,
        "Y_plus": 7.58,
        "Y_minus": 25.33,
    },
    # [REFERENCE: LNS, arXiv:1804.04706, Eq. (3.29)]
    "K0": {
        "pdg_pair": "311,-311",
        "plus": "K0 p",
        "minus": "K0bar p",
        "P": 11.93,
        "Y_plus": 9.08,
        "Y_minus": 19.09,
    },
}


# Evaluate the full unregulated moving eta factor used by raw mode
def raw_eta_factor(alpha_t: float, tau: int) -> complex:
    if tau not in {-1, 1}:
        raise ValueError(f"Invalid Regge signature tau={tau}")
    nearest_spin = round(alpha_t)
    spin_tau = 1 if abs(nearest_spin) % 2 == 0 else -1
    if alpha_t == nearest_spin and spin_tau == tau:
        return complex(math.inf, math.inf)
    phase = cmath.exp(-0.5j * math.pi * alpha_t)
    if tau == 1:
        return -phase / math.sin(0.5 * math.pi * alpha_t)
    return -1j * phase / math.cos(0.5 * math.pi * alpha_t)


# Evaluate the reduced rotating eta factor used by MRegge
def rotating_eta_factor(alpha_t: float, tau: int) -> complex:
    multiplier = {-1: -1j, 1: -1.0}.get(tau)
    if multiplier is None:
        raise ValueError(f"Invalid Regge signature tau={tau}")
    return multiplier * cmath.exp(-0.5j * math.pi * alpha_t)


# Validate one exact eta mode name
def _eta_mode(value: object, context: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{context} must be a string")
    if value not in {"raw", "rotating_t0", "rotating"}:
        raise ValueError(f"{context} has unknown mode {value}")
    return value


# Evaluate one configured MRegge eta factor
def eta_factor(alpha_t: float, alpha0: float, tau: int, eta_mode: str) -> complex:
    mode = _eta_mode(eta_mode, "Regge eta_mode")
    if mode == "raw":
        value = raw_eta_factor(alpha_t, tau)
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            return 0j
        return value
    if mode == "rotating_t0":
        return rotating_eta_factor(alpha0, tau)
    if mode == "rotating":
        return rotating_eta_factor(alpha_t, tau)
    raise ValueError(f"Regge eta_mode has unknown mode {mode}")


# Load the selected tune's mapped SOFT exchanges and signature modes
def load_config(path: Path = GENERAL_CARD, reader=None) -> dict[str, object]:
    general = read_card(path, reader)
    regge = general["PARAM_REGGE"]
    configuration = soft_exchange.load(general)
    soft_model = configuration.model_name
    soft_parameters = configuration.model
    soft_exchanges = soft_parameters["EXCHANGE"]
    exchange_def = configuration.definitions
    semantic_rows = {}
    for row in configuration.regge_exchanges:
        definition = exchange_def[row.soft_exchange]
        if row.role == "pomeron":
            key = "P"
        elif row.role == "odderon":
            key = "O"
        elif definition["crossing"] == 1:
            key = "C+"
        else:
            key = "C-"
        if key in semantic_rows:
            raise ValueError(f"PARAM_REGGE has multiple {key} rows for the DL table")
        semantic_rows[key] = row
    required = {exchange.key for exchange in CONFIGURED_EXCHANGES}
    if set(semantic_rows) != required:
        raise ValueError("PARAM_REGGE must map P, C-even, C-odd, and O rows for the DL table")
    mapped_names = tuple(
        semantic_rows[exchange.key].soft_exchange for exchange in CONFIGURED_EXCHANGES
    )
    residues = tuple(
        soft_exchange.proton_vertex(soft_parameters, name, soft_model) for name in mapped_names
    )
    exchanges = tuple(
        ConfiguredExchange(
            exchange=CONFIGURED_EXCHANGES[index],
            soft_exchange=mapped_names[index],
            trajectory_mode=exchange_def[mapped_names[index]]["trajectory_mode"],
            a0=soft_exchange.finite_number(
                soft_exchanges[mapped_names[index]]["alpha"][0],
                f"PARAM_SOFT.EXCHANGE.{mapped_names[index]}.alpha[0]",
            ),
            ap=soft_exchange.finite_number(
                soft_exchanges[mapped_names[index]]["alpha"][1],
                f"PARAM_SOFT.EXCHANGE.{mapped_names[index]}.alpha[1]",
            ),
            B=residues[index].cross_section_slope_per_gev2,
            beam_residue_t0=residues[index].beam_residue_per_gev,
            beam_residue_source=residues[index].source,
            tau=soft_exchange.exact_integer(
                exchange_def[mapped_names[index]]["tau"],
                f"PARAM_SOFT.EXCHANGE_DEF.{mapped_names[index]}.tau",
            ),
            eta_mode=_eta_mode(
                soft_exchanges[mapped_names[index]]["eta_mode"],
                f"PARAM_SOFT.EXCHANGE.{mapped_names[index]}.eta_mode",
            ),
        )
        for index in range(len(mapped_names))
    )
    for parameters in exchanges:
        eta_factor(
            parameters.a0,
            parameters.a0,
            parameters.tau,
            parameters.eta_mode,
        )
    photoprod_eta_mode = _eta_mode(regge["photoprod_eta_mode"], "PARAM_REGGE.photoprod_eta_mode")
    if photoprod_eta_mode not in {"rotating_t0", "rotating"}:
        raise ValueError("PARAM_REGGE.photoprod_eta_mode must be rotating_t0 or rotating")
    try:
        card_label = str(path.relative_to(ROOT))
    except ValueError:
        card_label = str(path)
    return {
        "card": card_label,
        "exchanges": exchanges,
        "soft_model": soft_model,
        "soft_eta_mode_P": _eta_mode(
            soft_exchanges[semantic_rows["P"].soft_exchange]["eta_mode"],
            f"PARAM_SOFT.EXCHANGE.{semantic_rows['P'].soft_exchange}.eta_mode",
        ),
        "soft_eta_mode_O": _eta_mode(
            soft_exchanges[semantic_rows["O"].soft_exchange]["eta_mode"],
            f"PARAM_SOFT.EXCHANGE.{semantic_rows['O'].soft_exchange}.eta_mode",
        ),
        "soft_eta_mode_3P": _eta_mode(soft_parameters["3P"]["eta_mode"], "PARAM_SOFT.3P.eta_mode"),
        "photoprod_eta_mode": photoprod_eta_mode,
    }


# Serialize one complex amplitude factor for JSON output
def complex_payload(value: complex) -> dict[str, float]:
    return {"real": value.real, "imag": value.imag}


# Serialize the mapped SOFT trajectories, proton couplings and signature factors
def exchange_output(
    configuration: dict[str, object] | None = None,
) -> list[dict[str, object]]:
    if configuration is None:
        configuration = load_config()
    rows = []
    for parameters in configuration["exchanges"]:
        eta = eta_factor(
            parameters.a0,
            parameters.a0,
            parameters.tau,
            parameters.eta_mode,
        )
        rows.append(
            {
                "exchange": parameters.exchange.key,
                "label": parameters.exchange.label,
                "soft_exchange": parameters.soft_exchange,
                "trajectory_mode": parameters.trajectory_mode,
                "a0": parameters.a0,
                "ap_GeV_minus2": parameters.ap,
                "B_GeV_minus2": parameters.B,
                "B_definition": "d ln |F_p(t)|^2 / dt at t=0 for the physical SOFT proton residue",
                "beam_residue_t0_GeV_minus1": parameters.beam_residue_t0,
                "beam_residue_source": parameters.beam_residue_source,
                "tau": parameters.tau,
                "eta_mode": parameters.eta_mode,
                "configured_eta_factor_t0": complex_payload(eta),
            }
        )
    return rows


# Split the original DL total-cross-section coefficients into exchange terms
def coefficients() -> dict[str, dict[str, float]]:
    out = {}
    for hadron, coeff in DL_MB.items():
        c_plus = 0.5 * (coeff["Y_plus"] + coeff["Y_minus"])
        c_minus = 0.5 * (coeff["Y_minus"] - coeff["Y_plus"])
        out[hadron] = {"P": coeff["P"], "C+": c_plus, "C-": c_minus}
    return out


# Compute the live physical-proton coupling for each fitted exchange
def beam_couplings(
    configuration: dict[str, object],
) -> dict[str, float]:
    out = {}
    for parameters in configuration["exchanges"]:
        if parameters.exchange not in EXCHANGES:
            continue
        residue = parameters.beam_residue_t0
        if not math.isfinite(residue) or residue <= 0.0:
            raise ValueError(
                f"Configured {parameters.exchange.key} beam residue must be finite and positive"
            )
        out[parameters.exchange.key] = residue
    required = {exchange.key for exchange in EXCHANGES}
    if set(out) != required:
        raise ValueError("Configured SOFT model is missing a DL beam residue")
    return out


# Compute the forward optical-theorem weight of each configured exchange
def optical_factors(
    configuration: dict[str, object],
) -> dict[str, float]:
    out = {}
    for parameters in configuration["exchanges"]:
        if parameters.exchange not in EXCHANGES:
            continue
        eta = eta_factor(
            parameters.a0,
            parameters.a0,
            parameters.tau,
            parameters.eta_mode,
        )
        weight = abs(eta.imag)
        if not math.isfinite(weight) or weight <= 0.0:
            raise ValueError(
                f"Configured {parameters.exchange.key} optical weight must be finite and positive"
            )
        out[parameters.exchange.key] = weight
    required = {exchange.key for exchange in EXCHANGES}
    if set(out) != required:
        raise ValueError("Configured SOFT model is missing a DL optical weight")
    return out


# Compute the symmetric reference couplings of the original DL convention
def reference_couplings(
    coeff_mb: dict[str, dict[str, float]],
) -> dict[tuple[str, str], float]:
    out = {}
    for exchange in EXCHANGES:
        exchange_key = exchange.key
        proton_coefficient = coeff_mb["p"][exchange_key]
        if not math.isfinite(proton_coefficient) or proton_coefficient <= 0.0:
            raise ValueError(f"DL {exchange_key} proton coefficient must be finite and positive")
        g_proton = math.sqrt(proton_coefficient * MB_TO_GEV2)
        for hadron in coeff_mb:
            coeff_gev2 = coeff_mb[hadron][exchange_key] * MB_TO_GEV2
            if not math.isfinite(coeff_gev2) or coeff_gev2 < 0.0:
                raise ValueError(
                    f"DL {exchange_key} {hadron} coefficient must be finite and nonnegative"
                )
            out[(hadron, exchange_key)] = coeff_gev2 / g_proton
    return out


# Normalize bare hadron couplings to the live physical SOFT proton coupling
def bare_couplings(
    coeff_mb: dict[str, dict[str, float]],
    beam_residues: dict[str, float],
) -> dict[tuple[str, str], float]:
    out = {}
    for exchange in EXCHANGES:
        key = exchange.key
        if key not in beam_residues:
            raise ValueError(f"Missing configured {key} beam residue")
        g_proton = beam_residues[key]
        if not math.isfinite(g_proton) or g_proton <= 0.0:
            raise ValueError(f"Configured {key} beam residue must be finite and positive")
        proton_coefficient = coeff_mb["p"][key]
        if not math.isfinite(proton_coefficient) or proton_coefficient <= 0.0:
            raise ValueError(f"DL {key} proton coefficient must be finite and positive")
        for hadron in coeff_mb:
            coefficient = coeff_mb[hadron][key]
            if not math.isfinite(coefficient) or coefficient < 0.0:
                raise ValueError(f"DL {key} {hadron} coefficient must be finite and nonnegative")
            out[(hadron, key)] = g_proton * coefficient / proton_coefficient
    return out


# Predict bare Born coefficients from the configured phase and coupling scale
def born_coefficients(
    residues: dict[tuple[str, str], float],
    optical_weights: dict[str, float],
) -> dict[str, dict[str, float]]:
    out = {hadron: {} for hadron, _ in residues}
    for exchange in EXCHANGES:
        key = exchange.key
        if key not in optical_weights:
            raise ValueError(f"Missing configured {key} optical weight")
        weight = optical_weights[key]
        if not math.isfinite(weight) or weight <= 0.0:
            raise ValueError(f"Configured {key} optical weight must be finite and positive")
        g_proton = residues[("p", key)]
        for hadron in out:
            out[hadron][key] = weight * g_proton * residues[(hadron, key)] / MB_TO_GEV2
    return out


# Compute effective non-strange and strange quark increments from mesons
def quark_couplings(
    residues: dict[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    out = {}
    for exchange in EXCHANGES:
        g_n = 0.5 * residues[("pi", exchange.key)]
        g_s = residues[("K", exchange.key)] - g_n
        out[("n", exchange.key)] = g_n
        out[("s", exchange.key)] = g_s
    return out


# Derive all continuum-row magnitudes using only DL and quark/isospin rules
def derived_couplings(
    residues: dict[tuple[str, str], float],
) -> dict[tuple[str, str], tuple[float, str]]:
    quark = quark_couplings(residues)
    out = {}

    for exchange in EXCHANGES:
        key = exchange.key
        g_phi = 2.0 * quark[("s", key)]
        for final_state in FINAL_STATES:
            if key == "C-" and final_state.key in {"pi0", "rho", "phi"}:
                out[(final_state.key, key)] = (0.0, "C parity zero")
            elif final_state.key == "n":
                out[(final_state.key, key)] = (
                    residues[("p", key)],
                    "isospin magnitude: n = p",
                )
            elif final_state.key == "pi0":
                out[(final_state.key, key)] = (
                    residues[("pi", key)],
                    "isospin: pi0 = pi",
                )
            elif final_state.key == "p":
                out[(final_state.key, key)] = (
                    residues[("p", key)],
                    "SOFT physical residue",
                )
            elif final_state.key == "K0":
                out[(final_state.key, key)] = (residues[("K0", key)], "isospin: K0 p = K+ n")
            elif final_state.direct_hadron is not None:
                out[(final_state.key, key)] = (
                    residues[(final_state.direct_hadron, key)],
                    "DL ratio to proton",
                )
            elif key == "C-" and final_state.key == "generic":
                out[(final_state.key, key)] = (
                    0.0,
                    "C-odd fallback set to zero",
                )
            elif final_state.key == "rho":
                out[(final_state.key, key)] = (
                    residues[("pi", key)],
                    "spin indep.: rho = pi",
                )
            elif final_state.key == "phi":
                out[(final_state.key, key)] = (
                    g_phi,
                    "additive ssbar",
                )
            elif final_state.key == "generic":
                out[(final_state.key, key)] = (
                    0.5 * (residues[("pi", key)] + g_phi),
                    "mean nnbar and ssbar residue",
                )
            else:
                raise ValueError(f"Unhandled final-state key {final_state.key}")
    return out


# Compute every model-specific PDG realization of one DL exchange
def model_exchange_pdgs(model: str, exchange: Exchange) -> tuple[int, ...]:
    return MODEL_EXCHANGE_PDGS[model][exchange.key]


# Compute whether one DL final state has an allowed lowest LS operator
def final_state_allowed(pdg: int, final_state: FinalState) -> bool:
    return not (pdg in ODD_SPIN_EXCHANGE_PDGS and final_state.key == "pi0")


# Build compact continuum card coupling block
def model_card_block(
    model: str,
    couplings: dict[tuple[str, str], tuple[float, str]],
    precision: int,
) -> dict[str, object]:
    exchange_rows = {}
    for exchange in EXCHANGES:
        for pdg in model_exchange_pdgs(model, exchange):
            rows = []
            for final_state in FINAL_STATES:
                if not final_state_allowed(pdg, final_state):
                    continue
                value = round(couplings[(final_state.key, exchange.key)][0], precision)
                rows.append([final_state.pdg[0], final_state.pdg[1], value, 0.0])
            exchange_rows[str(pdg)] = rows
    return {"g": exchange_rows}


# Compute the crossing-related DL coupling for one continuum particle family
def _family_coupling(rows: list[list[object]], pdg: list[int]) -> float:
    key = sorted(abs(value) for value in pdg)
    fallback = None
    for row in rows:
        if row[0] is None:
            fallback = float(row[2])
        elif sorted((abs(int(row[0])), abs(int(row[1])))) == key:
            return float(row[2])
    if fallback is None:
        raise ValueError(f"generated continuum couplings have no row for {pdg}")
    return fallback


# Compute the compact helicity row's parity and identical-particle orbit
def _helicity_orbit(row, parity: bool, identical: bool) -> set[tuple[float, float, int]]:
    h1, h2 = float(row[0]), float(row[1])
    m = soft_exchange.exact_integer(row[2], "Regge projection m") if len(row) == 5 else 0
    orbit = {(h1, h2, m)}
    if identical:
        orbit.add((h2, h1, m))
    if parity:
        orbit |= {(-a, -b, -c) for a, b, c in orbit}
    return orbit


# Compute the physical pole and allowed crossed LS tensors for one charge sector
def _pole_reference(exchange, pair, sector, block, particles):
    requested = (int(exchange), abs(int(pair[0])), abs(int(pair[1])) * (-1 if sector == "opposite" else 1))
    missing = [pdg for pdg in requested if pdg not in particles]
    if missing:
        raise ValueError(f"No embedded particle metadata for PDG {missing}; put explicit spinX2/P/C in GENERAL.PARAM_PDG")
    mother = particles[int(exchange)]
    first = particles[abs(int(pair[0]))]
    second_pdg = abs(int(pair[1])) * (-1 if sector == "opposite" else 1)
    second = particles[second_pdg]
    cp = block.get("CP")
    if not isinstance(cp, list) or len(cp) != 2 or any(type(flag) is not bool for flag in cp):
        raise ValueError("Continuum CP must contain [C,P] symmetry flags")
    allowed = sorted(spinmath.ls_rows(mother, first, second, *cp, crossed=True))
    # An empty operator space is a selection-rule zero, not a missing model
    # In particular, identical spin-zero pairs cannot couple to odd J
    return mother, first, second, allowed


# Convert a scalar DL/SOFT coupling to a declared physical-pole tensor norm
def _pole_norm_per_residue(pole, normalization: str) -> float:
    if normalization not in VERTEX_NORMALIZATIONS:
        raise ValueError(f"Unknown vertex normalization {normalization!r}")
    mother, first, second, allowed = pole
    if abs(first.pdg) != abs(second.pdg) or first.spin2 != second.spin2 or first.spin2 < 0:
        raise ValueError("DL pole-rms requires an elastic equal-spin hadron family")
    if first.pdg == 22 or second.pdg == 22:
        raise ValueError("DL pole-rms is for massive hadrons, not HERA photon vertices")
    return math.sqrt(first.spin2 + 1.0)


# Count physical helicities once using the C++ compact-row expansion
def _helicity_norm_weights(rows, anchor, block, model, identical, label, pole=None):
    parity = model == "GP" or block["CP"][1]
    orbits = [_helicity_orbit(rows[i], parity, identical) for i in anchor]
    seen, weights = {}, []
    if pole is not None:
        mother, first, second, _ = pole
        # GP pole spins have the signature parity; this is the same phase as
        # CrossedParityPhase / CrossedLegExchangePhase in the C++ parser
        spin_phase = (-1) ** ((first.spin2 + second.spin2 - mother.spin2) // 2)
        parity_phase = mother.parity * first.parity * second.parity * spin_phase
        exchange_phase = spin_phase
    for i, orbit in zip(anchor, orbits, strict=True):
        overlap = orbit & seen.keys()
        if overlap and model == "GP":
            raise ValueError(f"{label}.helicity repeats a symmetry orbit forbidden by C++")
        magnitude = soft_exchange.finite_number(rows[i][-2], f"{label}.magnitude")
        if pole is None:
            # Compatibility for callers needing only multiplicities
            if any(not math.isclose(magnitude, abs(seen[key]), rel_tol=1e-9, abs_tol=1e-9)
                   for key in overlap):
                raise ValueError(f"{label}.helicity has inconsistent symmetry-partner magnitudes")
            expanded = dict.fromkeys(orbit, complex(magnitude))
        else:
            a, b = float(rows[i][0]), float(rows[i][1])
            m = soft_exchange.exact_integer(rows[i][2], f"{label}.m") if model == "GP" else 0
            for h, part in ((a, first), (b, second)):
                if abs(h) > part.spin2 / 2 + 1e-9 or not spinmath.is_integer(h + part.spin2 / 2):
                    raise ValueError(f"{label}.helicity projection {h} is invalid for spin {part.spin2}/2")
            value = cmath.rect(magnitude, soft_exchange.finite_number(rows[i][-1], f"{label}.phase"))
            entries = [((a, b, m), value)]
            if identical:
                entries.append(((b, a, m), exchange_phase * value))
            if parity:
                entries += [((-x, -y, -z), parity_phase * v) for (x,y,z), v in list(entries)]
            expanded = {}
            for key, value in entries:
                if key in expanded and abs(expanded[key] - value) > 1e-9 * max(1., abs(value)):
                    raise ValueError(f"{label}.helicity has a nonzero symmetry-forbidden component {key}")
                expanded[key] = value
            for key in overlap:
                if abs(expanded[key] - seen[key]) > 1e-9 * max(1., abs(expanded[key]), abs(seen[key])):
                    raise ValueError(f"{label}.helicity has inconsistent symmetry-partner phases or magnitudes")
        weights.append(len(orbit - seen.keys()))
        seen.update(expanded)
    return orbits, weights


# Zero sector magnitudes, retaining phases
def _zero_coupling_updates(base, block, model):
    label = f"CON_{model}.json:" + ".".join(base)
    basis = block.get("basis")
    if basis == "crossed_auto_min_L" and model in {"MP", "XP"}:
        g = block.get("g")
        if not isinstance(g, list) or len(g) != 2 or any(not math.isfinite(float(x)) for x in g):
            raise ValueError(f"{label}.g must contain finite [magnitude,phase]")
        return [ScalarUpdate((*base, "g", 0), f"{label}.g.magnitude", 0.0)]
    if basis not in {"crossed_ls", "crossed_helicity"}:
        raise ValueError(f"{label} has an unsupported basis {basis!r}")
    field = "g_ls" if basis == "crossed_ls" else "helicity"
    rows, size = block.get(field), 5 if model == "GP" else 4
    if not isinstance(rows, list) or not rows or any(not isinstance(r, list) or len(r) != size for r in rows):
        raise ValueError(f"{label}.{field} must have {size}-column rows")
    return [ScalarUpdate((*base, field, i, size-2), f"{label}.{field}[{i}].magnitude", 0.0)
            for i in range(len(rows))]


# Reduced pole tensor norm, using m=0 for GP
def pole_tensor_norm(block, pole, model):
    mother, first, second, allowed = pole
    if not allowed:
        return 0.0
    basis = block["basis"]
    if basis == "crossed_auto_min_L":
        return abs(float(block["g"][0])) * spinmath.ls_norm(mother, first, second, *allowed[0])
    rows = block["g_ls" if basis == "crossed_ls" else "helicity"]
    anchor = [i for i, row in enumerate(rows) if model != "GP" or row[2] == 0]
    if basis == "crossed_ls":
        weights = [spinmath.ls_norm(mother, first, second, *rows[i][:2]) ** 2 for i in anchor]
    else:
        identical = first.pdg == second.pdg and spinmath.has_cparity(first) and spinmath.has_cparity(second)
        _, weights = _helicity_norm_weights(rows, anchor, block, model, identical, "pole norm", pole=pole)
    return math.sqrt(math.fsum(w * float(rows[i][-2]) ** 2 for i, w in zip(anchor, weights, strict=True)))


# Normalize every representation to the same declared pole strength
def _coupling_updates(base, block, target, pole, model, shape, *, normalization=DEFAULT_VERTEX_NORMALIZATION):
    if not math.isfinite(target) or target < 0.0:
        raise ValueError("DL residue must be finite and nonnegative")
    if normalization not in VERTEX_NORMALIZATIONS:
        raise ValueError(f"Unknown vertex normalization {normalization!r}")
    if shape not in {"preserve", "lowest"}:
        raise ValueError(f"Unknown coupling shape {shape!r}")
    if target == 0.0 or (pole is not None and not pole[3]):
        return _zero_coupling_updates(base, block, model)
    basis = block.get("basis")
    label = f"CON_{model}.json:" + ".".join(base)
    if basis == "crossed_auto_min_L" and model in {"MP", "XP"}:
        g = block.get("g")
        if not isinstance(g, list) or len(g) != 2 or any(not math.isfinite(float(x)) for x in g) or g[0] < 0:
            raise ValueError(f"{label}.g must contain finite [nonnegative magnitude,phase]")
        mother, first, second, allowed = pole
        raw = spinmath.ls_norm(mother, first, second, *allowed[0])
        value = target * _pole_norm_per_residue(pole, normalization) / raw
        updates = [ScalarUpdate((*base, "g", 0), f"{label}.g.magnitude", round(value, COUPLING_PRECISION))]
        if shape == "lowest":
            updates.append(ScalarUpdate((*base, "g", 1), f"{label}.g.phase", 0.0))
        return updates
    if basis not in {"crossed_ls", "crossed_helicity"}:
        raise ValueError(f"{label} requires crossed_ls, crossed_helicity or crossed_auto_min_L")
    field = "g_ls" if basis == "crossed_ls" else "helicity"
    rows = block.get(field)
    size = 5 if model == "GP" else 4
    if not isinstance(rows, list) or not rows or any(not isinstance(row, list) or len(row) != size for row in rows):
        raise ValueError(f"{label}.{field} must be a nonempty array of {size}-column rows")
    for row in rows:
        if any(not math.isfinite(soft_exchange.finite_number(value, label)) for value in row) or float(row[-2]) < 0.0:
            raise ValueError(f"{label}.{field} has an invalid coupling row")
    anchor = [i for i, row in enumerate(rows) if size == 4 or soft_exchange.exact_integer(row[2], "m") == 0]
    if not anchor:
        raise ValueError(f"{label}.{field} requires an m=0 pole anchor")
    mother, first, second, allowed = pole
    lowest = allowed[0]
    pole_norm = spinmath.ls_norm(mother, first, second, *lowest)
    target_norm = target * _pole_norm_per_residue(pole, normalization)
    seed_coefficient = target_norm / pole_norm
    if field == "g_ls":
        seen = set()
        for row in rows:
            L = soft_exchange.exact_integer(row[0], f"{label}.L")
            S = soft_exchange.exact_integer(row[1], f"{label}.S")
            if (L, S) not in allowed:
                raise ValueError(f"{label}.g_ls contains forbidden pole LS tensor {(L,S)}")
            key = (L, S, abs(soft_exchange.exact_integer(row[2], "m")) if model == "GP" else 0)
            if key in seen:
                raise ValueError(f"{label}.g_ls repeats an LS or m-reflection orbit")
            seen.add(key)
        weights = [spinmath.ls_norm(mother, first, second, *rows[i][:2]) ** 2 for i in anchor]
    else:
        identical = first.pdg == second.pdg and spinmath.has_cparity(first) and spinmath.has_cparity(second)
        # Validate all supplied GP m columns, including their parity partners
        if model == "GP":
            _helicity_norm_weights(rows, list(range(len(rows))), block, model, identical, label, pole=pole)
        orbits, weights = _helicity_norm_weights(rows, anchor, block, model, identical, label, pole=pole)
    norm = math.sqrt(math.fsum(weight * float(rows[i][-2]) ** 2
                               for i, weight in zip(anchor, weights, strict=True)))
    if shape == "preserve" and norm == 0.0:
        # There is no nonzero shape to preserve. Initialize the allowed lowest
        # operator using the existing layout instead of requesting more input
        shape = "lowest"
    if shape == "preserve":
        rounding = 0.5 * 10.0 ** -COUPLING_PRECISION * math.sqrt(math.fsum(weights))
        scale = 1.0 if abs(norm-target_norm) <= rounding else target_norm / norm
        initial = [float(row[-2]) * scale for row in rows]
    elif field == "g_ls":
        if not any(tuple(rows[i][:2]) == lowest for i in anchor):
            raise ValueError(f"{label}.g_ls omits the lowest allowed LS tensor {lowest} at m=0")
        initial = [seed_coefficient if i in anchor and tuple(row[:2]) == lowest else 0.0
                   for i,row in enumerate(rows)]
    else:
        helicity = spinmath.helicity_rows(mother, first, second, allowed, lowest, p_symmetry=False)
        values = {(row[0], row[1], 0): seed_coefficient * cmath.rect(row[2], row[3]) for row in helicity}
        missing = {key for key,value in values.items() if abs(value)>1e-12} - set().union(*orbits)
        if missing:
            raise ValueError(f"{label}.helicity omits lowest-LS helicities: {sorted(missing)}")
        initial = [values.get((row[0],row[1],0),0.0) if i in anchor else 0.0 for i,row in enumerate(rows)]
    updates = []
    for i, (row,value) in enumerate(zip(rows, initial, strict=True)):
        coordinate = ",".join(str(item) for item in row[:-2])
        row_label = f"{label}.{field}[{coordinate}]"
        updates.append(ScalarUpdate((*base, field, i, size-2), f"{row_label}.magnitude",
                                    round(abs(value), COUPLING_PRECISION)))
        if shape == "lowest":
            phase = -math.pi if complex(value).real < 0.0 else 0.0
            updates.append(ScalarUpdate((*base, field, i, size-1), f"{row_label}.phase", phase))
    return updates


# Build scalar updates for every allowed crossing-related charge sector
def card_updates(
    card: dict[str, object],
    generated: dict[str, object],
    model: str,
    tune_dir: Path | None = None,
    *,
    shape: str,
    normalization: str = DEFAULT_VERTEX_NORMALIZATION,
    reader=None,
) -> list[ScalarUpdate]:
    if normalization not in VERTEX_NORMALIZATIONS:
        raise ValueError(f"Unknown vertex normalization {normalization!r}")
    if shape not in {"preserve", "lowest"}:
        raise ValueError("DL coupling shape must be preserve or lowest")
    source_dir = ROOT / "modeldata" / "TUNE0" if tune_dir is None else tune_dir
    particles = load_particles(source_dir, reader)
    if model == "GP":
        _, regge = load_regge(source_dir, reader)
        for group, spin in zip(regge.groups, regge.pole_spin, strict=True):
            for pdg in group:
                if pdg in particles and particles[pdg].spin2 < 0:
                    particles[pdg] = replace(particles[pdg], spin2=2 * spin)
    updates = []
    for exchange, generated_rows in generated["g"].items():
        if exchange not in card:
            continue
        exchange_block = card[exchange]
        if not isinstance(exchange_block, dict):
            raise ValueError(f"CON_{model}.json exchange {exchange} must be an object")
        for pair_key, pair_block in exchange_block.items():
            pair = json.loads(pair_key)
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError(f"CON_{model}.json has invalid pair key {pair_key}")
            sectors = [name for name in ("same", "opposite", "self") if name in pair_block]
            unsupported = set(pair_block) - {*sectors, "FF_transfer", "FF_offshell", "reggeize", "pveto"}
            if unsupported or not sectors:
                raise ValueError(f"CON_{model}.json:{exchange}.{pair_key} has invalid sectors")
            if "self" in sectors and len(sectors) != 1:
                raise ValueError(
                    f"CON_{model}.json:{exchange}.{pair_key} self sector cannot coexist "
                    "with same or opposite sectors"
                )
            target = _family_coupling(generated_rows, pair)
            for sector_name in sectors:
                block = pair_block[sector_name]
                pole = _pole_reference(exchange, pair, sector_name, block, particles)
                updates.extend(_coupling_updates(
                    (exchange, pair_key, sector_name), block, target, pole, model, shape,
                    normalization=normalization,
                ))
    if not updates:
        raise ValueError(f"CON_{model}.json contains no supported DL target channels")
    return updates


# Check rounded coupling norms for each sector
def normalization_audit(card, revised, generated, model, tune_dir, normalization, reader=None):
    particles = load_particles(tune_dir, reader)
    records = []
    for exchange, generated_rows in generated["g"].items():
        if exchange not in card:
            continue
        for pair_key, family in card[exchange].items():
            pair = json.loads(pair_key)
            scalar = _family_coupling(generated_rows, pair)
            for sector in ("same", "opposite", "self"):
                if sector not in family:
                    continue
                old = family[sector]
                new = revised[exchange][pair_key][sector]
                pole = _pole_reference(exchange, pair, sector, old, particles)
                mother, first, second, allowed = pole
                forbidden = not allowed
                target = 0.0 if forbidden else scalar
                field = "g" if new["basis"] == "crossed_auto_min_L" else (
                    "g_ls" if new["basis"] == "crossed_ls" else "helicity")
                magnitudes = [new["g"][0]] if field == "g" else [r[-2] for r in new[field]]
                if target == 0:
                    if any(x != 0 for x in magnitudes):
                        raise ValueError(f"{model}:{exchange}:{pair_key}:{sector}: nonzero forbidden coupling")
                    actual, expected, bound = 0.0, 0.0, 0.0
                else:
                    actual = pole_tensor_norm(new, pole, model)
                    expected = target * _pole_norm_per_residue(pole, normalization)
                    if field == "g":
                        sensitivity = spinmath.ls_norm(mother, first, second, *allowed[0])
                    elif field == "g_ls":
                        sensitivity = math.sqrt(math.fsum(
                            spinmath.ls_norm(mother, first, second, *r[:2]) ** 2
                            for r in new[field] if model != "GP" or r[2] == 0))
                    else:
                        anchor = [i for i,r in enumerate(new[field]) if model != "GP" or r[2] == 0]
                        identical = first.pdg == second.pdg and spinmath.has_cparity(first) and spinmath.has_cparity(second)
                        _, weights = _helicity_norm_weights(new[field], anchor, new, model, identical,
                                                           "normalization audit", pole=pole)
                        sensitivity = math.sqrt(math.fsum(weights))
                    bound = 0.50001 * 10.0 ** -COUPLING_PRECISION * sensitivity + 1e-12 * max(1.0, expected)
                    if abs(actual-expected) > bound:
                        raise ValueError(f"{model}:{exchange}:{pair_key}:{sector}: rounded pole norm "
                                         f"{actual} disagrees with target {expected}")
                records.append({
                    "model":model, "exchange":exchange, "pair":pair, "sector":sector,
                    "basis":new["basis"], "J_pole":mother.spin2/2,
                    "incoming_spin_states":first.spin2+1,
                    "scalar_residue":target,
                    "target_pole_norm":expected, "actual_pole_norm":actual,
                    "absolute_error":abs(actual-expected), "rounding_bound":bound,
                    "status":"selection_rule_zero" if forbidden else ("zero_residue" if target==0 else "matched"),
                })
    return records


# Preview and push selected DL continuum couplings
def push_cards(
    general_path: Path,
    card_blocks: dict[str, object],
    *,
    confirm=None,
    shape: str = DEFAULT_SHAPE,
    normalization: str = DEFAULT_VERTEX_NORMALIZATION,
) -> bool:
    targets = {model: general_path.parent / f"CON_{model}.json" for model in card_blocks}
    reader = CardReader({general_path, *targets.values()})
    updates = {path: card_updates(
        reader.read(path), card_blocks[model], model, tune_dir=general_path.parent,
        shape=shape, normalization=normalization, reader=reader) for model, path in targets.items()}
    return push_json5_updates(updates, confirm=confirm, reader=reader)


# Build a machine-readable payload containing inputs, derivations and card rows
def output(
    coeff_mb: dict[str, dict[str, float]],
    reference_residues: dict[tuple[str, str], float],
    bare_residues: dict[tuple[str, str], float],
    quark: dict[tuple[str, str], float],
    couplings: dict[tuple[str, str], tuple[float, str]],
    models: list[str],
    precision: int,
    configuration: dict[str, object] | None = None,
    *,
    normalization: str = DEFAULT_VERTEX_NORMALIZATION,
) -> dict[str, object]:
    if configuration is None:
        configuration = load_config()
    configured_rows = exchange_output(configuration)
    beam_residues = beam_couplings(configuration)
    optical_weights = optical_factors(configuration)
    bare_coefficients = born_coefficients(bare_residues, optical_weights)
    return {
        "convention": {
            "sigma_hplus_p": "P*s^eps + (C+ - C-)*s^-eta",
            "sigma_hminus_p": "P*s^eps + (C+ + C-)*s^-eta",
            "epsilon": DL_EPSILON,
            "eta": DL_ETA,
            "mb_to_GeV_minus2": MB_TO_GEV2,
        },
        "tool_version": TOOL_VERSION,
        "spin_normalization": {
            "observable": "unpolarized total cross section: linear forward-amplitude trace",
            "optical_average": "sum(Im M_forward_diagonal) / ((2*s_h+1)*2*s)",
            "scalar_residue_average": "cancels for spin-independent helicity-conserving residues",
            "extra_scalar_spin_factor": 1.0,
            "vertex_normalization": normalization,
            "pole_seed": "||H_pole||^2/(2*s_h+1) = g_h^2",
            "normalized_shape": "H_pole = g_h*sqrt(2*s_h+1)*T/||T||_F",
            "model_closure": "spin-averaged reduced-pole strength; distinct from the optical forward trace",
            "scope": "additional pole seed convention, not TP matching; GP anchor m=0 at pole J",
            "CON_report": "actual proposed target cards; full scalar targets are in residue_seed_rows",
            "massive_vector_spin_average": 3,
            "fermion_spin_average": 2,
            "reggeon_spin_average": "none",
            "shape_default": DEFAULT_SHAPE,
            "zero_shape": "initialize lowest allowed LS in the existing layout",
            "forbidden_sector": "zero existing magnitudes; retain layout",
            "scalar_scale": "physical SOFT proton residue, DL hadron/proton ratios",
            "vector_model": "spin-independent rho=pi; phi=2*K-pi; all three massive helicities",
            "nonzero_GP_m": "common rescaling with m=0; unchanged relative complex coefficients",
        },
        "configured_exchange_parameters": configured_rows,
        "eta_convention": {
            "general_card": configuration["card"],
            "configured_eta_mode": [row["eta_mode"] for row in configured_rows],
            "raw_formula": "-(1+tau*exp(-i*pi*alpha))/sin(pi*alpha)",
            "rotating_even_formula": "-exp(-i*pi*alpha/2)",
            "rotating_odd_formula": "-i*exp(-i*pi*alpha/2)",
            "reference_pomeron": "i",
            "reference_even": "-cot(pi*alpha/2) + i",
            "reference_odd": "-tan(pi*alpha/2) - i",
            "soft_modes": {
                "model": configuration["soft_model"],
                "eta_mode_P": configuration["soft_eta_mode_P"],
                "eta_mode_O": configuration["soft_eta_mode_O"],
                "3P": {"eta_mode": configuration["soft_eta_mode_3P"]},
                "scope": (
                    'Independent Pomeron and Odderon modes, 3P.eta_mode for high mass SD/DD'
                ),
            },
            "photoprod_eta_mode": configuration["photoprod_eta_mode"],
        },
        "coefficients_mb": coeff_mb,
        "coefficients_GeV_minus2": {
            hadron: {key: coefficients[key] * MB_TO_GEV2 for key in ("P", "C+", "C-")}
            for hadron, coefficients in coeff_mb.items()
        },
        "residue_factorization": {
            "formula": ("g_p^E = g_p^E(SOFT,t=0), g_h^E = g_p^E*C_hp^E/C_pp^E"),
            "beam_residues_GeV_minus1": beam_residues,
            "optical_weights": optical_weights,
            "bare_Born_coefficients_mb": bare_coefficients,
            "bare_to_DL_coefficient_scale": {
                exchange.key: bare_coefficients["p"][exchange.key] / coeff_mb["p"][exchange.key]
                for exchange in EXCHANGES
            },
            "scope": (
                'SOFT fixes g_p, DL fixes g_h/g_p'
            ),
            "secondary_assumption": (
                "the effective C-even f2/a2 and C-odd rho/omega sums are each assigned "
                "to one mapped SOFT exchange"
            ),
        },
        "residue_couplings_GeV_minus1": [
            {
                "exchange": exchange.key,
                "hadron": hadron,
                "symmetric_DL_reference": reference_residues[(hadron, exchange.key)],
                "bare_value": bare_residues[(hadron, exchange.key)],
            }
            for exchange in EXCHANGES
            for hadron in coeff_mb
        ],
        "effective_quark_couplings_GeV_minus1": [
            {
                "exchange": exchange.key,
                "g_n": quark[("n", exchange.key)],
                "g_s": quark[("s", exchange.key)],
            }
            for exchange in EXCHANGES
        ],
        "continuum_couplings_GeV_minus1": [
            {
                "exchange": exchange.key,
                "final_state": final_state.label,
                "pdg": list(final_state.pdg),
                "value": couplings[(final_state.key, exchange.key)][0],
                "origin": couplings[(final_state.key, exchange.key)][1],
            }
            for exchange in EXCHANGES
            for final_state in FINAL_STATES
        ],
        "CON": {model: model_card_block(model, couplings, precision) for model in models},
    }


# Print the coefficient, coupling and derived continuum tables
def print_tables(
    coeff_mb: dict[str, dict[str, float]],
    reference_residues: dict[tuple[str, str], float],
    bare_residues: dict[tuple[str, str], float],
    quark: dict[tuple[str, str], float],
    couplings: dict[tuple[str, str], tuple[float, str]],
    configuration: dict[str, object] | None = None,
    *,
    normalization: str = DEFAULT_VERTEX_NORMALIZATION,
) -> None:
    print("Donnachie-Landshoff continuum couplings")
    print("  sigma(h+ p) = P*s^eps + (C+ - C-)*s^-eta")
    print("  sigma(h- p) = P*s^eps + (C+ + C-)*s^-eta")
    print(f"  eps = {DL_EPSILON:.4f}, eta = {DL_ETA:.4f}, s in GeV^2")
    print(f"  1 mb = {MB_TO_GEV2:.5f} GeV^-2\n")

    if configuration is None:
        configuration = load_config()
    _print_exchanges(configuration)
    _print_hadrons(coeff_mb, reference_residues, bare_residues, configuration)
    _print_quarks(quark, bare_residues)
    _print_continuum(couplings, bare_residues, normalization)


# Print SOFT trajectories and proton couplings
def _print_exchanges(configuration):
    configured_rows = []
    for parameters in configuration["exchanges"]:
        eta = eta_factor(
            parameters.a0,
            parameters.a0,
            parameters.tau,
            parameters.eta_mode,
        )
        configured_rows.append(
            [
                parameters.exchange.key,
                parameters.soft_exchange,
                parameters.trajectory_mode,
                f"{parameters.a0:.4f}",
                f"{parameters.ap:.4f}",
                f"{parameters.B:.4f}",
                f"{parameters.beam_residue_t0:.4f}",
                f"{parameters.tau:+d}",
                parameters.eta_mode,
                f"{eta.real:.4f} {eta.imag:+.4f}i",
            ]
        )
    print(f"SOFT model: {configuration['soft_model']}")
    print(f"Input card: {configuration['card']}")
    common.print_table(
        "Trajectories and proton couplings at t=0",
        [
            "key",
            "SOFT exchange",
            "trajectory",
            "alpha0",
            "alpha' [GeV^-2]",
            "B_p(0) [GeV^-2]",
            "g_p(0) [GeV^-1]",
            "tau",
            "eta_mode",
            "eta(0)",
        ],
        configured_rows,
        right_align={3, 4, 5, 6, 7, 9},
    )
    print("  B_p = d ln|F_p|^2/dt at t=0, including the Good Walker projection")
    print("  C_hp,bare = |Im eta(0)| g_h g_p [GeV^-2]")


# Print DL coefficients and hadron couplings
def _print_hadrons(coeff_mb, reference_residues, bare_residues, configuration):
    optical_weights = optical_factors(configuration)
    bare_coefficients = born_coefficients(bare_residues, optical_weights)
    coefficient_rows = []
    for hadron, coeff in coeff_mb.items():
        coefficient_rows.append(
            [hadron, DL_MB[hadron]["pdg_pair"], *(f"{coeff[key]:.4f}" for key in ("P", "C+", "C-"))]
        )
    common.print_table(
        "DL coefficients [mb]", ["hadron", "PDG pair", "P", "C+", "C-"], coefficient_rows, right_align={2, 3, 4}
    )

    residue_rows = []
    for exchange in EXCHANGES:
        for hadron in coeff_mb:
            residue_rows.append(
                [
                    exchange.label,
                    exchange.key,
                    hadron,
                    f"{coeff_mb[hadron][exchange.key]:.4f}",
                    f"{coeff_mb[hadron][exchange.key] * MB_TO_GEV2:.4f}",
                    f"{reference_residues[(hadron, exchange.key)]:.4f}",
                    f"{bare_residues[(hadron, exchange.key)]:.4f}",
                    f"{bare_coefficients[hadron][exchange.key]:.4f}",
                ]
            )
    common.print_table(
        "\nHadron couplings from SOFT and DL",
        [
            "exchange",
            "key",
            "hadron",
            "C_hp [mb]",
            "C_hp [GeV^-2]",
            "g_h(DL reference) [GeV^-1]",
            "g_h(SOFT bare) [GeV^-1]",
            "C_hp,bare [mb]",
        ],
        residue_rows,
        right_align={3, 4, 5, 6, 7},
        group_by=0,
    )
    print("  DL reference: g_p=sqrt(C_pp), |Im eta|=1")
    print("  SOFT fixes bare g_p, DL fixes g_h/g_p=C_hp/C_pp")
    print("  Bare Born coefficients need not equal effective DL coefficients")


# Print quark couplings
def _print_quarks(quark, bare_residues):
    quark_rows = []
    for exchange in EXCHANGES:
        quark_rows.append(
            [
                exchange.label,
                exchange.key,
                f"{quark[('n', exchange.key)]:.4f}",
                f"{quark[('s', exchange.key)]:.4f}",
                f"{3.0 * quark[('n', exchange.key)]:.4f}",
                f"{bare_residues[('p', exchange.key)]:.4f}",
            ]
        )
    common.print_table(
        "\nQuark couplings [GeV^-1]: g_pi=2g_n, g_K=g_n+g_s",
        ["exchange", "key", "g_n", "g_s", "3g_n", "g_p(SOFT bare)"],
        quark_rows,
        right_align={2, 3, 4, 5},
    )


# Print continuum couplings and normalization
def _print_continuum(couplings, bare_residues, normalization):
    pomeron_ppbar = bare_residues[("p", "P")]
    continuum_rows = []
    for exchange in EXCHANGES:
        for final_state in FINAL_STATES:
            value, origin = couplings[(final_state.key, exchange.key)]
            continuum_rows.append(
                [
                    exchange.label,
                    exchange.key,
                    final_state.label,
                    "*" if final_state.pdg[0] is None else f"{final_state.pdg[0]},{final_state.pdg[1]}",
                    origin,
                    f"{value:.4f}",
                    f"{value / pomeron_ppbar:.4f}",
                ]
            )
    common.print_table(
        "\nContinuum couplings",
        ["exchange", "key", "final state", "PDG pair", "origin", "g [GeV^-1]", "g/g_Ppp(SOFT bare)"],
        continuum_rows,
        right_align={5, 6},
        group_by=0,
    )
    print(f"\nContinuum vertex normalization: {normalization}")
    print("  Pole convention: ||H||^2/(2*s_h+1)=g_h^2, GP normalized at m=0")
    print("  Additional pole assumption beyond the optical theorem")
    print("  Preserve mode rescales all components, lowest mode initializes lowest LS")
    print("  Zero anchors use lowest LS, forbidden sectors stay zero")
    print("  Each C parity Reggeon sum seeds one effective SOFT exchange")


# Parse the common output and continuum-model options
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Derive continuum couplings from DL fits')
    parser.add_argument(
        "--format", choices=("table", "json", "cards"), default="table",
        help='table with update preview, or JSON/cards output',
    )
    parser.add_argument(
        "--model",
        choices=("all", "MP", "XP", "GP"),
        default="all",
        help="continuum model(s) included in card output",
    )
    parser.add_argument(
        "--precision", type=int, default=6, help="decimal places used in generated card couplings"
    )
    parser.add_argument(
        "--tune-general",
        type=Path,
        default=GENERAL_CARD,
        help="GENERAL.json providing SOFT inputs and the directory of the continuum cards",
    )
    parser.add_argument(
        "--push",
        action="store_true",
        help="show the tables and confirm continuum card updates, regardless of --format",
    )
    parser.add_argument(
        "--shape", choices=("preserve", "lowest"), default=DEFAULT_SHAPE,
        help="preserve relative amplitudes and phases by default, or initialize the lowest LS tensor",
    )
    parser.add_argument("--dry-run", action="store_true", help="preview target coefficient changes without writing")
    parser.add_argument("--yes", action="store_true", help="confirm writes, only together with --push")
    return parser.parse_args()


# Compute couplings and update the selected continuum cards
def main() -> int:
    args = parse_args()
    if args.precision < 1 or args.precision > 15:
        raise ValueError("--precision must be between 1 and 15")
    if args.yes and not args.push:
        raise ValueError("--yes requires --push")
    general_path = args.tune_general.resolve()
    if general_path.name != "GENERAL.json":
        raise ValueError("--tune-general must name GENERAL.json")

    models = ["MP", "XP", "GP"] if args.model == "all" else [args.model]
    targets = {model: general_path.parent / f"CON_{model}.json" for model in models}
    reader = CardReader({general_path, *targets.values()})
    configuration = load_config(general_path, reader)
    coeff_mb = coefficients()
    reference_residues = reference_couplings(coeff_mb)
    bare_residues = bare_couplings(coeff_mb, beam_couplings(configuration))
    quark = quark_couplings(bare_residues)
    couplings = derived_couplings(bare_residues)
    payload = output(coeff_mb, reference_residues, bare_residues, quark, couplings,
                            models, args.precision, configuration, normalization=DEFAULT_VERTEX_NORMALIZATION)

    shape = args.shape
    updates, proposed, checks = {}, {}, []

    for model, path in targets.items():
        card = reader.read(path)

        # Keep internal calculations unrounded; --precision only formats scalar
        # summary rows. Written coupling magnitudes retain nine decimal places
        seeds = model_card_block(model, couplings, 15)
        changes = card_updates(card, seeds, model, tune_dir=general_path.parent,
                                          shape=shape, normalization=DEFAULT_VERTEX_NORMALIZATION, reader=reader)
        updates[path] = changes
        revised = copy.deepcopy(card)
        for change in changes:
            _set(revised, change.path, change.value)
        proposed[model] = revised
        checks.extend(normalization_audit(card, revised, seeds, model, general_path.parent,
                                          DEFAULT_VERTEX_NORMALIZATION, reader))

    # Preflight every reference/conflict even when no write is requested
    _prepare_updates(updates, reader)
    payload["residue_seed_rows"] = payload.pop("CON")
    payload["CON"] = proposed
    payload["spin_normalization"]["shape"] = shape
    payload["inputs"] = [str(general_path), *(str(path) for path in targets.values())]
    payload["normalization_checks"] = checks
    if args.push or args.format == "table" or args.dry_run:
        print_tables(coeff_mb, reference_residues, bare_residues, quark, couplings,
                     configuration, normalization=DEFAULT_VERTEX_NORMALIZATION)
        print(f"\nSelected coupling shape: {shape}")
        print(f"Normalization audit: {len(checks)} sectors checked before writing")
        for model in models:
            selected = [r for r in checks if r["model"] == model]
            zero = sum(r["status"] == "selection_rule_zero" for r in selected)
            error = max((r["absolute_error"] for r in selected), default=0.0)
            print(f"  {model}: {len(selected)} sectors, {zero} selection-rule zeros, max norm error {error:.3g}")
        print(f"\nTarget-card preview: {general_path.parent}")
        push_json5_updates(updates, confirm=(lambda rows: True) if args.yes else None, dry_run=args.dry_run, reader=reader)
    elif args.format == "json":
        print(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        print(common.dumps({"spin_normalization": payload["spin_normalization"], "CON": proposed}))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError, OSError) as exc:
        print(f"DL_couplings.py: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
