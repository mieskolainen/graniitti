# Numba accelerated FastJet sequential recombination jet clustering
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import numba
import numpy as np


# Compute transverse momentum squared from Cartesian components
@numba.njit(cache=True)
def jet_pt2(px, py):
    return px * px + py * py


# Compute transverse momentum from Cartesian components
@numba.njit(cache=True)
def jet_pt(px, py):
    return math.hypot(px, py)


# Compute pseudorapidity with a finite value for nonzero transverse momentum
@numba.njit(cache=True)
def jet_eta(px, py, pz):
    pt = jet_pt(px, py)
    if pt > 0.0:
        return math.asinh(pz / pt)
    if abs(pz) > 0.0:
        return math.copysign(math.inf, pz)
    return 0.0


# Compute longitudinal rapidity with stable forward numerical behaviour
@numba.njit(cache=True)
def jet_rapidity(px, py, pz, energy):
    pt2 = jet_pt2(px, py)
    transverse_mass2 = (energy - abs(pz)) * (energy + abs(pz))
    transverse_mass2 = max(transverse_mass2, pt2)
    if transverse_mass2 > 0.0:
        return math.asinh(pz / math.sqrt(transverse_mass2))
    if abs(pz) > 0.0:
        return math.copysign(math.inf, pz)
    return 0.0


# Compute azimuth in the interval from minus pi to pi
@numba.njit(cache=True)
def jet_phi(px, py):
    return math.atan2(py, px)


# Compute invariant mass with negative numerical mass squared mapped to zero
@numba.njit(cache=True)
def jet_mass(px, py, pz, energy):
    mass2 = energy * energy - px * px - py * py - pz * pz
    return math.sqrt(max(mass2, 0.0))


# Compute wrapped absolute azimuthal separation
@numba.njit(cache=True)
def abs_delta_phi(first, second):
    delta = (first - second + math.pi) % (2.0 * math.pi) - math.pi
    return abs(delta)


# Compute anti-kT angular distance squared
@numba.njit(cache=True)
def angular_distance2(rapidity_first, phi_first, rapidity_second, phi_second):
    delta_rapidity = rapidity_first - rapidity_second
    delta_phi = abs_delta_phi(phi_first, phi_second)
    return delta_rapidity * delta_rapidity + delta_phi * delta_phi


# Compute the flattened tile index for one rapidity and azimuth coordinate
@numba.njit(cache=True)
def tiled_index(rapidity, phi, rapidity_min, rapidity_tile_size, phi_tile_size, rapidity_tile_count, phi_tile_count):
    rapidity_tile = int(math.floor((rapidity - rapidity_min) / rapidity_tile_size))
    rapidity_tile = min(max(rapidity_tile, 0), rapidity_tile_count - 1)
    phi_tile = int(math.floor((phi + math.pi) / phi_tile_size)) % phi_tile_count
    return rapidity_tile * phi_tile_count + phi_tile


# Insert one active pseudojet at the head of its spatial tile
@numba.njit(cache=True)
def tile_insert(index, tile, tile_head, node_tile, tile_previous, tile_next):
    previous_head = tile_head[tile]
    node_tile[index] = tile
    tile_previous[index] = -1
    tile_next[index] = previous_head
    if previous_head >= 0:
        tile_previous[previous_head] = index
    tile_head[tile] = index


# Remove one active pseudojet from its spatial tile
@numba.njit(cache=True)
def tile_remove(index, tile_head, node_tile, tile_previous, tile_next):
    tile = node_tile[index]
    previous = tile_previous[index]
    following = tile_next[index]
    if previous >= 0:
        tile_next[previous] = following
    else:
        tile_head[tile] = following
    if following >= 0:
        tile_previous[following] = previous
    tile_previous[index] = -1
    tile_next[index] = -1


# Compute the nearest active pseudojet inside neighbouring spatial tiles
@numba.njit(cache=True)
def nearest_tiled_neighbor(index, node_rapidity, node_phi, node_tile, tile_head, tile_next, phi_tile_count, radius2):
    center_tile = node_tile[index]
    center_rapidity_tile = center_tile // phi_tile_count
    center_phi_tile = center_tile % phi_tile_count
    rapidity_tile_count = tile_head.size // phi_tile_count
    nearest = -1
    nearest_distance2 = radius2

    for rapidity_offset in range(-1, 2):
        rapidity_tile = center_rapidity_tile + rapidity_offset
        if rapidity_tile < 0 or rapidity_tile >= rapidity_tile_count:
            continue
        for phi_offset in range(-1, 2):
            phi_tile = (center_phi_tile + phi_offset) % phi_tile_count
            tile = rapidity_tile * phi_tile_count + phi_tile
            candidate = tile_head[tile]
            while candidate >= 0:
                if candidate != index:
                    distance2 = angular_distance2(
                        node_rapidity[index], node_phi[index], node_rapidity[candidate], node_phi[candidate]
                    )
                    if distance2 < nearest_distance2:
                        nearest = candidate
                        nearest_distance2 = distance2
                candidate = tile_next[candidate]
    return nearest, nearest_distance2


# Add the unique neighbouring tiles around one tile to a compact union
@numba.njit(cache=True)
def add_neighbor_tiles(center_tile, tile_union, union_count, rapidity_tile_count, phi_tile_count):
    center_rapidity_tile = center_tile // phi_tile_count
    center_phi_tile = center_tile % phi_tile_count
    for rapidity_offset in range(-1, 2):
        rapidity_tile = center_rapidity_tile + rapidity_offset
        if rapidity_tile < 0 or rapidity_tile >= rapidity_tile_count:
            continue
        for phi_offset in range(-1, 2):
            phi_tile = (center_phi_tile + phi_offset) % phi_tile_count
            tile = rapidity_tile * phi_tile_count + phi_tile
            duplicate = False
            for index in range(union_count):
                if tile_union[index] == tile:
                    duplicate = True
                    break
            if not duplicate:
                tile_union[union_count] = tile
                union_count += 1
    return union_count


# Update local nearest-neighbour information after removal or recombination
@numba.njit(cache=True)
def update_tiled_neighbors(
    removed_first,
    removed_second,
    merged,
    first_tile,
    second_tile,
    merged_tile,
    node_rapidity,
    node_phi,
    node_tile,
    nearest,
    nearest_distance2,
    tile_head,
    tile_next,
    rapidity_tile_count,
    phi_tile_count,
    tile_union,
    radius2,
):
    union_count = add_neighbor_tiles(first_tile, tile_union, 0, rapidity_tile_count, phi_tile_count)
    if removed_second >= 0:
        union_count = add_neighbor_tiles(second_tile, tile_union, union_count, rapidity_tile_count, phi_tile_count)
        union_count = add_neighbor_tiles(merged_tile, tile_union, union_count, rapidity_tile_count, phi_tile_count)

    for union_index in range(union_count):
        candidate = tile_head[tile_union[union_index]]
        while candidate >= 0:
            if nearest[candidate] == removed_first or (removed_second >= 0 and nearest[candidate] == removed_second):
                neighbor, distance2 = nearest_tiled_neighbor(
                    candidate, node_rapidity, node_phi, node_tile, tile_head, tile_next, phi_tile_count, radius2
                )
                nearest[candidate] = neighbor
                nearest_distance2[candidate] = distance2

            if merged >= 0 and candidate != merged:
                distance2 = angular_distance2(
                    node_rapidity[candidate], node_phi[candidate], node_rapidity[merged], node_phi[merged]
                )
                if distance2 < nearest_distance2[candidate]:
                    nearest[candidate] = merged
                    nearest_distance2[candidate] = distance2
                if distance2 < nearest_distance2[merged]:
                    nearest[merged] = candidate
                    nearest_distance2[merged] = distance2
            candidate = tile_next[candidate]


# Compute the anti-kT distance represented by one nearest-neighbour entry
@numba.njit(cache=True)
def antikt_candidate_distance(index, nearest, nearest_distance2, inverse_pt2, inverse_radius2):
    momentum_factor = inverse_pt2[index]
    neighbor = nearest[index]
    if neighbor >= 0:
        momentum_factor = min(momentum_factor, inverse_pt2[neighbor])
    return momentum_factor * (nearest_distance2[index] * inverse_radius2)


# Cluster E-scheme anti-kT jets with N2Tiled nearest-neighbour updates
# [REFERENCE: M. Cacciari, G. P. Salam and G. Soyez, JHEP 04 (2008) 063, arXiv:0802.1189]
# [REFERENCE: M. Cacciari, G. P. Salam and G. Soyez, Eur. Phys. J. C 72 (2012) 1896, arXiv:1111.6097]
@numba.njit(cache=True)
def cluster_antikt(px, py, pz, energy, radius, pt_min=0.0, abs_eta_max=1.0e6):
    if px.ndim != 1 or py.ndim != 1 or pz.ndim != 1 or energy.ndim != 1:
        raise ValueError("Four-momentum components must be one-dimensional arrays")
    particle_count = px.size
    if py.size != particle_count or pz.size != particle_count or energy.size != particle_count:
        raise ValueError("Four-momentum component arrays must have equal lengths")
    empty_jets = np.empty((0, 4), dtype=np.float64)
    empty_labels = np.empty(0, dtype=np.int64)
    radius2 = radius * radius
    if (
        particle_count == 0
        or not math.isfinite(radius)
        or radius <= 0.0
        or not math.isfinite(radius2)
        or radius2 < np.finfo(np.float64).tiny
        or math.isnan(pt_min)
        or pt_min < 0.0
        or math.isnan(abs_eta_max)
        or abs_eta_max <= 0.0
    ):
        return empty_jets, empty_labels
    inverse_radius2 = 1.0 / radius2

    maximum_nodes = 2 * particle_count
    node_p4 = np.zeros((maximum_nodes, 4), dtype=np.float64)
    node_inverse_pt2 = np.zeros(maximum_nodes, dtype=np.float64)
    node_rapidity = np.zeros(maximum_nodes, dtype=np.float64)
    node_phi = np.zeros(maximum_nodes, dtype=np.float64)
    active = np.zeros(maximum_nodes, dtype=np.bool_)
    head = np.full(maximum_nodes, -1, dtype=np.int64)
    tail = np.full(maximum_nodes, -1, dtype=np.int64)
    next_particle = np.full(particle_count, -1, dtype=np.int64)
    roots = np.full(particle_count, -1, dtype=np.int64)

    active_count = 0
    rapidity_min = math.inf
    rapidity_max = -math.inf
    for index in range(particle_count):
        if not (
            math.isfinite(px[index])
            and math.isfinite(py[index])
            and math.isfinite(pz[index])
            and math.isfinite(energy[index])
        ):
            continue
        if energy[index] <= 0.0:
            continue
        rapidity = jet_rapidity(px[index], py[index], pz[index], energy[index])
        # Omit particles along the beam with infinite rapidity
        if not math.isfinite(rapidity):
            continue
        node_p4[index, 0] = px[index]
        node_p4[index, 1] = py[index]
        node_p4[index, 2] = pz[index]
        node_p4[index, 3] = energy[index]
        # Bound the inverse for vanishing transverse momentum without dropping constituents
        node_inverse_pt2[index] = 1.0 / max(jet_pt2(px[index], py[index]), 1.0e-300)
        node_rapidity[index] = rapidity
        node_phi[index] = jet_phi(px[index], py[index])
        active[index] = True
        head[index] = index
        tail[index] = index
        active_count += 1
        rapidity_min = min(rapidity_min, node_rapidity[index])
        rapidity_max = max(rapidity_max, node_rapidity[index])

    next_node = particle_count
    root_count = 0
    if active_count == 0:
        return empty_jets, np.full(particle_count, -1, dtype=np.int64)

    rapidity_tile_size = max(0.1, radius)
    phi_tile_count = max(3, int(math.floor(2.0 * math.pi / rapidity_tile_size)))
    phi_tile_size = 2.0 * math.pi / phi_tile_count
    minimum_rapidity_tile = int(math.floor(rapidity_min / rapidity_tile_size))
    maximum_rapidity_tile = int(math.floor(rapidity_max / rapidity_tile_size))
    rapidity_tile_count = maximum_rapidity_tile - minimum_rapidity_tile + 1
    tiled_rapidity_min = minimum_rapidity_tile * rapidity_tile_size
    tile_head = np.full(rapidity_tile_count * phi_tile_count, -1, dtype=np.int64)
    node_tile = np.full(maximum_nodes, -1, dtype=np.int64)
    tile_previous = np.full(maximum_nodes, -1, dtype=np.int64)
    tile_next = np.full(maximum_nodes, -1, dtype=np.int64)
    nearest = np.full(maximum_nodes, -1, dtype=np.int64)
    nearest_distance2 = np.full(maximum_nodes, radius2, dtype=np.float64)
    tile_union = np.empty(27, dtype=np.int64)

    for index in range(particle_count):
        if not active[index]:
            continue
        tile = tiled_index(
            node_rapidity[index],
            node_phi[index],
            tiled_rapidity_min,
            rapidity_tile_size,
            phi_tile_size,
            rapidity_tile_count,
            phi_tile_count,
        )
        tile_insert(index, tile, tile_head, node_tile, tile_previous, tile_next)

    for index in range(particle_count):
        if not active[index]:
            continue
        neighbor, distance2 = nearest_tiled_neighbor(
            index, node_rapidity, node_phi, node_tile, tile_head, tile_next, phi_tile_count, radius2
        )
        nearest[index] = neighbor
        nearest_distance2[index] = distance2

    while active_count > 0:
        minimum_distance = math.inf
        minimum_first = -1
        for index in range(next_node):
            if not active[index]:
                continue
            neighbor = nearest[index]
            if neighbor >= 0 and not active[neighbor]:
                neighbor, distance2 = nearest_tiled_neighbor(
                    index, node_rapidity, node_phi, node_tile, tile_head, tile_next, phi_tile_count, radius2
                )
                nearest[index] = neighbor
                nearest_distance2[index] = distance2
            distance = antikt_candidate_distance(index, nearest, nearest_distance2, node_inverse_pt2, inverse_radius2)
            if distance < minimum_distance:
                minimum_distance = distance
                minimum_first = index

        if minimum_first < 0:
            break
        minimum_second = nearest[minimum_first]
        if minimum_second < 0:
            first_tile = node_tile[minimum_first]
            tile_remove(minimum_first, tile_head, node_tile, tile_previous, tile_next)
            active[minimum_first] = False
            active_count -= 1
            roots[root_count] = minimum_first
            root_count += 1
            update_tiled_neighbors(
                minimum_first,
                -1,
                -1,
                first_tile,
                -1,
                -1,
                node_rapidity,
                node_phi,
                node_tile,
                nearest,
                nearest_distance2,
                tile_head,
                tile_next,
                rapidity_tile_count,
                phi_tile_count,
                tile_union,
                radius2,
            )
            continue

        merged = next_node
        next_node += 1
        node_p4[merged] = node_p4[minimum_first] + node_p4[minimum_second]
        merged_pt2 = jet_pt2(node_p4[merged, 0], node_p4[merged, 1])
        node_inverse_pt2[merged] = 1.0 / max(merged_pt2, 1.0e-300)
        node_rapidity[merged] = jet_rapidity(
            node_p4[merged, 0], node_p4[merged, 1], node_p4[merged, 2], node_p4[merged, 3]
        )
        node_phi[merged] = jet_phi(node_p4[merged, 0], node_p4[merged, 1])
        head[merged] = head[minimum_first]
        tail[merged] = tail[minimum_second]
        next_particle[tail[minimum_first]] = head[minimum_second]

        first_tile = node_tile[minimum_first]
        second_tile = node_tile[minimum_second]
        tile_remove(minimum_first, tile_head, node_tile, tile_previous, tile_next)
        tile_remove(minimum_second, tile_head, node_tile, tile_previous, tile_next)
        active[minimum_first] = False
        active[minimum_second] = False
        active[merged] = True
        active_count -= 1

        merged_tile = tiled_index(
            node_rapidity[merged],
            node_phi[merged],
            tiled_rapidity_min,
            rapidity_tile_size,
            phi_tile_size,
            rapidity_tile_count,
            phi_tile_count,
        )
        tile_insert(merged, merged_tile, tile_head, node_tile, tile_previous, tile_next)
        update_tiled_neighbors(
            minimum_first,
            minimum_second,
            merged,
            first_tile,
            second_tile,
            merged_tile,
            node_rapidity,
            node_phi,
            node_tile,
            nearest,
            nearest_distance2,
            tile_head,
            tile_next,
            rapidity_tile_count,
            phi_tile_count,
            tile_union,
            radius2,
        )

    root_pt2 = np.empty(root_count, dtype=np.float64)
    for index in range(root_count):
        root_pt2[index] = jet_pt2(node_p4[roots[index], 0], node_p4[roots[index], 1])
    root_order = np.argsort(-root_pt2)
    sorted_roots = np.empty(root_count, dtype=np.int64)
    for index in range(root_count):
        sorted_roots[index] = roots[root_order[index]]

    selected_count = 0
    for index in range(root_count):
        root = sorted_roots[index]
        if jet_pt(node_p4[root, 0], node_p4[root, 1]) < pt_min:
            continue
        if abs(jet_eta(node_p4[root, 0], node_p4[root, 1], node_p4[root, 2])) > abs_eta_max:
            continue
        selected_count += 1

    jets = np.empty((selected_count, 4), dtype=np.float64)
    labels = np.full(particle_count, -1, dtype=np.int64)
    selected_index = 0
    for index in range(root_count):
        root = sorted_roots[index]
        if jet_pt(node_p4[root, 0], node_p4[root, 1]) < pt_min:
            continue
        if abs(jet_eta(node_p4[root, 0], node_p4[root, 1], node_p4[root, 2])) > abs_eta_max:
            continue
        jets[selected_index] = node_p4[root]
        constituent = head[root]
        while constituent >= 0:
            labels[constituent] = selected_index
            constituent = next_particle[constituent]
        selected_index += 1

    return jets, labels
