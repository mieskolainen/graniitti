# Data and MC covariance handling
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import os

import numpy as np
from core.io.files import ensure_dir
from core.io.serialize import load_json_file
from core.io.serialize import sha256_file as file_sha256
from core.numerics import array
from core.stats import hist, uncertainty
from core.stats.uncertainty import source_covariance, stable_sum, stable_sum_squares
from scipy import sparse

DATA_COVARIANCE_ARRAYS = (
    "data_stat_uncertainty",
    "data_count_stat_uncertainty",
    "data_block_values",
    "data_block_columns",
    "data_block_row_offsets",
    "data_block_shape",
    "data_collective_vectors",
    "data_component_labels",
) + tuple(f"{name}_{field}" for name in ("stat_local", "syst_local", "mc_cross")
          for field in ("values", "columns", "row_offsets"))

DATA_COVARIANCE_SCHEMA = 1
STRUCTURED_COVARIANCE_KEYS = (
    "data_block_values",
    "data_block_columns",
    "data_block_row_offsets",
    "data_block_shape",
    "data_collective_vectors",
    "data_component_labels",
)


# Materialize complete covariance matrices only when explicitly requested for small comparisons
class CovariancePayload(dict):
    # Derive dense diagnostic matrices without storing them in initialization or worker state
    def __missing__(self, key):
        if key not in {"data_stat_correlation", "data_stat_covariance", "data_syst_covariance",
                       "data_total_covariance", "mc_cross_correlation"}:
            raise KeyError(key)
        return covariance_selection(self, key)


# Encode a SciPy sparse covariance using portable numerical arrays
def sparse_covariance_arrays(matrix, name):
    matrix = sparse.csr_matrix(matrix)
    matrix.eliminate_zeros()
    return {f"{name}_values": matrix.data.copy(), f"{name}_columns": matrix.indices.astype(np.int64),
            f"{name}_row_offsets": matrix.indptr.astype(np.int64)}


# Restore one sparse covariance with the common retained-bin shape
def sparse_covariance_matrix(payload, name):
    shape = tuple(np.asarray(payload["data_block_shape"], dtype=np.int64))
    return sparse.csr_matrix((payload[f"{name}_values"], payload[f"{name}_columns"],
                              payload[f"{name}_row_offsets"]), shape=shape)


# Compute only the requested covariance rows and columns including collective sources
def covariance_selection(payload, name="data_total_covariance", indices=None):
    if name in payload:
        matrix = np.asarray(payload[name], dtype=float)
        return matrix if indices is None else matrix[np.ix_(indices, indices)]
    size = int(payload["data_block_shape"][0])
    indices = np.arange(size) if indices is None else np.asarray(indices, dtype=np.int64)
    if name == "mc_cross_correlation":
        return sparse_covariance_matrix(payload, "mc_cross")[indices][:, indices].toarray()
    structure = prepare_covariance_structure(payload)
    vectors = np.asarray(structure["collective_vectors"])[indices]
    if name == "data_total_covariance":
        block = structure["matrix"][indices][:, indices].toarray()
    else:
        statistical = name in {"data_stat_covariance", "data_stat_correlation"}
        if name not in {"data_stat_covariance", "data_stat_correlation", "data_syst_covariance"}:
            raise KeyError(name)
        block = sparse_covariance_matrix(payload, "stat_local" if statistical else "syst_local")[indices][:, indices].toarray()
        selected = np.asarray([source["category"] == "statistical" for source in payload["data_collective_sources"]], dtype=bool)
        vectors = vectors[:, selected if statistical else ~selected]
        if statistical:
            sigma = np.asarray(payload["data_count_stat_uncertainty"])[indices]
            cross = sparse_covariance_matrix(payload, "mc_cross")[indices][:, indices].toarray()
            block += float(payload["mc_correlation_scale"]) * cross * np.outer(sigma, sigma)
    block += vectors @ vectors.T
    if name == "data_stat_correlation":
        return uncertainty.covariance_to_correlation(block, np.asarray(payload["data_stat_uncertainty"])[indices])
    return block


# Compute deterministic valid-bin layout records for nested tune histograms
def build_layout(data: list) -> list[dict]:
    layout = []
    offset = 0
    for dataset_index, dataset in enumerate(data):
        if dataset is None:
            continue
        for subset_index, subset in enumerate(dataset):
            for observable in sorted(subset):
                histogram = subset[observable]["hdata"]
                valid = np.asarray(histogram.valid, dtype=bool)
                bins = np.flatnonzero(valid)
                sources = subset[observable].get("uncertainties", [])
                mc_correlation_model = (
                    "normalized"
                    if any(
                        source["category"] == "statistical" and bool(source.get("mc_correlation_eligible", False))
                        for source in sources
                    )
                    else "counting"
                )
                layout.append(
                    {
                        "dataset": dataset_index,
                        "subset": subset_index,
                        "observable": observable,
                        "mc_correlation_model": mc_correlation_model,
                        "bins": bins.tolist(),
                        "start": offset,
                        "stop": offset + len(bins),
                    }
                )
                offset += len(bins)
    return layout


# Compute the global vector indices represented by one layout record
def record_indices(record: dict) -> np.ndarray:
    return np.arange(record["start"], record["stop"], dtype=np.int64)


# Copy one local source covariance into a global covariance matrix
def add_local_covariance(target: np.ndarray, record: dict, covariance: np.ndarray) -> None:
    bins = np.asarray(record["bins"], dtype=np.int64)
    indices = record_indices(record)
    target[np.ix_(indices, indices)] += covariance[np.ix_(bins, bins)]


# Assemble local covariance terms and retain collective nuisance vectors
def _data_covariance_terms(data: list, layout: list[dict], *, sparse_output: bool = False) -> dict:
    size = 0 if not layout else layout[-1]["stop"]
    local_blocks = {"statistical": [], "systematic": []}
    data_count_stat_variance = np.zeros(size, dtype=float)
    collective: dict[tuple[str, str], np.ndarray] = {}
    source_metadata = []

    for record in layout:
        item = data[record["dataset"]][record["subset"]][record["observable"]]
        sources = item.get("uncertainties", [])
        bins = np.asarray(record["bins"], dtype=np.int64)
        indices = record_indices(record)
        blocks = {category: sparse.csr_matrix((len(bins), len(bins))) for category in local_blocks}
        for source in sources:
            source_metadata.append(
                {
                    "dataset": record["dataset"],
                    "subset": record["subset"],
                    "observable": record["observable"],
                    "name": source["name"],
                    "category": source["category"],
                    "correlation": source["correlation"],
                    "effect": source.get("effect", "additive"),
                    "scope": source.get("scope"),
                    "provenance": source.get("provenance"),
                }
            )
            if source["correlation"] == "collective":
                local_covariance = None
            elif source["correlation"] == "uncorrelated":
                sigma = 0.5 * (np.abs(np.asarray(source["up"])) + np.abs(np.asarray(source["down"])))
                local_covariance = sparse.diags(sigma[bins] ** 2, format="csr")
            else:
                local_covariance = sparse.csr_matrix(source_covariance(source)[np.ix_(bins, bins)])
            if source["category"] == "statistical" and (
                source["correlation"] == "uncorrelated" or bool(source.get("mc_correlation_eligible", False))
            ):
                variance = (local_covariance.diagonal() if local_covariance is not None else
                            np.diag(source_covariance(source))[bins])
                data_count_stat_variance[indices] += np.clip(variance, 0.0, None)
            if source["correlation"] != "collective":
                category = "statistical" if source["category"] == "statistical" else "systematic"
                blocks[category] += local_covariance
                continue

            scope = source.get("scope")
            if not scope:
                raise ValueError(f"Collective uncertainty '{source['name']}' requires a scope")
            vector = np.asarray(
                source.get(
                    "shift",
                    0.5
                    * (np.abs(np.asarray(source["up"], dtype=float)) + np.abs(np.asarray(source["down"], dtype=float))),
                ),
                dtype=float,
            )
            global_vector = collective.setdefault((source["category"], scope), np.zeros(size, dtype=float))
            global_vector[indices] = vector[bins]

        for category, block in blocks.items():
            local_blocks[category].append(block)

    collective_sources = [{"category": category, "scope": scope} for category, scope in collective]
    collective_vectors = np.column_stack(tuple(collective.values())) if collective else np.empty((size, 0), dtype=float)
    local = {category: sparse.block_diag(blocks, format="csr") if blocks else sparse.csr_matrix((size, size))
             for category, blocks in local_blocks.items()}
    return {
        "data_stat_local": local["statistical"] if sparse_output else local["statistical"].toarray(),
        "data_syst_local": local["systematic"] if sparse_output else local["systematic"].toarray(),
        "data_count_stat_uncertainty": np.sqrt(data_count_stat_variance),
        "data_uncertainty_sources": source_metadata,
        "data_collective_sources": collective_sources,
        "data_collective_vectors": collective_vectors,
    }


# Form the covariance carried by selected collective nuisance vectors
def collective_covariance(vectors: np.ndarray, sources: list[dict], *, category: str | None = None) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != len(sources):
        raise ValueError("Collective covariance sources and vectors are incompatible")
    selected = np.ones(len(sources), dtype=bool)
    if category is not None:
        selected = np.asarray(
            [
                source.get("category") == "statistical"
                if category == "statistical"
                else source.get("category") != "statistical"
                for source in sources
            ],
            dtype=bool,
        )
    retained = vectors[:, selected]
    return retained @ retained.T


# Assemble complete data statistical and systematic covariance terms
def assemble_data_covariances(data: list, layout: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    terms = _data_covariance_terms(data, layout)
    vectors = terms["data_collective_vectors"]
    sources = terms["data_collective_sources"]
    data_stat_covariance = terms["data_stat_local"] + collective_covariance(vectors, sources, category="statistical")
    data_syst_covariance = terms["data_syst_local"] + collective_covariance(vectors, sources, category="systematic")
    return (
        data_stat_covariance,
        data_syst_covariance,
        terms["data_count_stat_uncertainty"],
        terms["data_uncertainty_sources"],
    )


# Build accepted-event bin assignments for every histogram in one subset
def subset_assignments(mcdata: dict, obs: dict) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    output = {}
    event_ids = np.asarray(mcdata["event_ids"], dtype=np.int64)
    for observable, config in obs.items():
        values = np.asarray(mcdata["data"][observable])
        bins = config["bins"]
        if not isinstance(bins, dict):
            output[observable] = (event_ids, hist.bin_indices(values, bins))
            continue

        members = hist.assign_2d_hyperbins(x=values[:, 0:2], hyperbins=bins, symmetrized=config["symmetrized_fill"])
        for subcategory in bins:
            x0_bins = bins[subcategory][0]
            x1_bins = bins[subcategory][1]
            x2_bins = bins[subcategory][2]
            name = f"{observable}_{hist.bins2txt(x0_bins)}_{hist.bins2txt(x1_bins)}"
            local = hist.bin_indices(values[:, 2], x2_bins)
            member_mask = np.zeros(len(local), dtype=bool)
            member_mask[members[subcategory]] = True
            local[~member_mask] = -1
            output[name] = (event_ids, local)
    return output


# Compute independently generated sample groups for one dataset
def mcdata_sample_groups(entry: object) -> list[dict]:
    if isinstance(entry, dict) and "sample_groups" in entry:
        return list(entry["sample_groups"])
    return [{"mcdata": entry, "sample": None, "set_indices": None}]


# Build a sparse event-by-bin indicator matrix for one generated dataset sample
def indicator_matrix(
    *, dataset_index: int, mcdata: list[dict], obs: list[dict], layout: list[dict], set_indices: list[int] | None = None
) -> tuple[sparse.csr_matrix, np.ndarray]:
    if set_indices is None:
        set_indices = list(range(len(mcdata)))
    if len(set_indices) != len(mcdata):
        raise ValueError("MC correlation set routing does not match the generated sample")
    local_by_set = {set_index: local for local, set_index in enumerate(set_indices)}
    records = [item for item in layout if item["dataset"] == dataset_index and item["subset"] in local_by_set]
    if not records:
        return sparse.csr_matrix((0, 0), dtype=float), np.empty(0, dtype=float)

    sample_ids = np.asarray(mcdata[0]["sample_event_ids"], dtype=np.int64)
    sample_weights = np.asarray(mcdata[0]["sample_weights"], dtype=float)
    if len(sample_ids) != len(sample_weights):
        raise ValueError("MC correlation event ids and weights have different lengths")
    if len(np.unique(sample_ids)) != len(sample_ids):
        raise ValueError("MC correlation event ids are not unique")
    if not np.all(np.isfinite(sample_weights)):
        raise ValueError("MC correlation event weights contain non-finite values")
    id_to_row = {int(event_id): index for index, event_id in enumerate(sample_ids)}
    column_map = {}
    for record in records:
        for position, bin_index in enumerate(record["bins"]):
            column_map[(record["subset"], record["observable"], int(bin_index))] = record["start"] + position

    rows = []
    columns = []
    assignments = {
        set_index: subset_assignments(mcdata[local_index], obs[set_index])
        for set_index, local_index in local_by_set.items()
    }
    for record in records:
        event_ids, local_bins = assignments[record["subset"]][record["observable"]]
        for event_id, local_bin in zip(event_ids, local_bins, strict=True):
            key = (record["subset"], record["observable"], int(local_bin))
            if local_bin < 0 or key not in column_map:
                continue
            rows.append(id_to_row[int(event_id)])
            columns.append(column_map[key])

    size = 0 if not layout else layout[-1]["stop"]
    values = np.ones(len(rows), dtype=float)
    matrix = sparse.coo_matrix((values, (rows, columns)), shape=(len(sample_ids), size)).tocsr()
    return matrix, sample_weights


# Estimate cross-observable correlations from one MC event sample
def estimate_mc_cross_correlation(
    *, obs: list, mcdata_by_dataset: list, layout: list[dict], sparse_output: bool = False
) -> np.ndarray | sparse.csr_matrix:
    size = 0 if not layout else layout[-1]["stop"]
    covariance = sparse.csr_matrix((size, size), dtype=float)
    for dataset_index, dataset_mcdata in enumerate(mcdata_by_dataset):
        if dataset_mcdata is None:
            continue
        for group in mcdata_sample_groups(dataset_mcdata):
            mcdata = group["mcdata"]
            if not mcdata:
                continue
            indicator, weights = indicator_matrix(
                dataset_index=dataset_index,
                mcdata=mcdata,
                obs=obs[dataset_index],
                layout=layout,
                set_indices=group["set_indices"],
            )
            if indicator.shape[0] == 0:
                continue
            weight_scale = float(np.max(np.abs(weights)))
            if not np.isfinite(weight_scale) or weight_scale <= 0.0:
                continue
            weights = weights / weight_scale
            weight_sum = stable_sum(weights)
            weight2_sum = stable_sum_squares(weights)
            if weight_sum <= 0.0 or weight2_sum <= 0.0:
                continue
            probability = np.asarray(indicator.T @ weights).ravel() / weight_sum
            centered = np.zeros(size, dtype=float)
            for record in layout:
                if record["dataset"] == dataset_index and record["mc_correlation_model"] == "normalized":
                    centered[record_indices(record)] = 1.0
            centered_probability = centered * probability
            weighted = indicator.multiply(weights[:, None])
            joint = (weighted.T @ weighted).tocsr()
            marginal2 = np.asarray(indicator.T @ (weights**2)).ravel()
            centered_column = sparse.csr_matrix(centered_probability[:, None])
            marginal_column = sparse.csr_matrix(marginal2[:, None])
            block = (joint - marginal_column @ centered_column.T - centered_column @ marginal_column.T
                     + weight2_sum * (centered_column @ centered_column.T))
            covariance += block

    diagonal = np.sqrt(np.clip(covariance.diagonal(), 0.0, None))
    inverse = sparse.diags(np.divide(1.0, diagonal, out=np.zeros_like(diagonal), where=diagonal > 0.0))
    correlation = inverse @ covariance @ inverse
    correlation = 0.5 * (correlation + correlation.T)
    if np.any(np.abs(correlation.data) > 1.0 + 1.0e-10):
        raise ValueError("MC event overlap produced an invalid correlation")
    correlation.data = np.clip(correlation.data, -1.0, 1.0)
    observable = np.empty(size, dtype=np.int64)
    for index, record in enumerate(layout):
        observable[record_indices(record)] = index
    entries = correlation.tocoo()
    retained = observable[entries.row] != observable[entries.col]
    correlation = sparse.csr_matrix((entries.data[retained], (entries.row[retained], entries.col[retained])),
                                    shape=(size, size))
    correlation.eliminate_zeros()
    return correlation if sparse_output else correlation.toarray()


# Summarize the MC samples used only for data statistical correlations
def mc_correlation_sample_metadata(
    mcdata_by_dataset: list, mc_correlation_datacards: list[dict] | None = None
) -> list[dict]:
    records = []
    for dataset_index, dataset_mcdata in enumerate(mcdata_by_dataset):
        if dataset_mcdata is None:
            continue
        card = mc_correlation_datacards[dataset_index] if mc_correlation_datacards is not None else {}
        for group in mcdata_sample_groups(dataset_mcdata):
            mcdata = group["mcdata"]
            if not mcdata:
                continue
            event_ids = np.asarray(mcdata[0]["sample_event_ids"], dtype=np.int64)
            weights = np.asarray(mcdata[0]["sample_weights"], dtype=float)
            if len(event_ids) != len(weights):
                raise ValueError("MC correlation event ids and weights have different lengths")
            weight_sum = stable_sum(weights)
            weight2_sum = stable_sum_squares(weights)
            weight_scale = float(np.max(np.abs(weights))) if len(weights) > 0 else 0.0
            scaled_weights = weights / weight_scale if weight_scale > 0.0 else weights
            scaled_sum = stable_sum(scaled_weights)
            scaled_sum2 = stable_sum_squares(scaled_weights)
            record = {
                "dataset": dataset_index,
                "event_count_requested": (int(card["nevents"]) if "nevents" in card else None),
                "event_count_read": int(len(event_ids)),
                "generator_weighted": (bool(card["weighted"]) if "weighted" in card else None),
                "event_weights_nonuniform": bool(
                    len(weights) > 0 and not np.allclose(weights, weights[0], rtol=1.0e-12, atol=0.0)
                ),
                "weight_sum": weight_sum,
                "weight_squared_sum": weight2_sum,
                "effective_event_count": (scaled_sum**2 / scaled_sum2 if scaled_sum2 > 0.0 else 0.0),
            }
            if group["sample"] is not None:
                record["sample"] = group["sample"]
            records.append(record)
    return records


# Compute the lowest covariance eigenvalue globally or over exact blocks
def minimum_covariance_eigenvalue(covariance: np.ndarray, component_labels: np.ndarray | None = None) -> float:
    covariance = np.asarray(covariance, dtype=float)
    if covariance.size == 0:
        return np.inf
    if component_labels is None:
        return float(np.min(np.linalg.eigvalsh(covariance)))
    labels = np.asarray(component_labels, dtype=np.int64)
    if labels.shape != (len(covariance),):
        raise ValueError("Covariance component labels have incompatible dimensions")
    minimum = np.inf
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        block = covariance[np.ix_(indices, indices)]
        minimum = min(minimum, float(np.min(np.linalg.eigvalsh(block))))
    return minimum


# Scale only MC-derived cross terms until the covariance is positive semidefinite
def scale_mc_cross_covariance_to_psd(
    data_stat_covariance: np.ndarray, mc_cross_covariance: np.ndarray, *, component_labels: np.ndarray | None = None,
    reference_scale: float | None = None,
) -> tuple[float, np.ndarray]:
    data_stat = np.asarray(data_stat_covariance, dtype=float)
    mc_cross = np.asarray(mc_cross_covariance, dtype=float)
    if data_stat.size == 0:
        return 1.0, data_stat.copy()
    scale = max(float(np.max(np.abs(np.diag(data_stat)))), 1.0e-30) if reference_scale is None else reference_scale
    tolerance = 1.0e-10 * scale
    if minimum_covariance_eigenvalue(data_stat, component_labels) < -tolerance:
        raise ValueError("Data statistical covariance is not positive semidefinite")

    candidate = 0.5 * (data_stat + mc_cross + (data_stat + mc_cross).T)
    if minimum_covariance_eigenvalue(candidate, component_labels) >= -tolerance:
        return 1.0, candidate

    lower = 0.0
    upper = 1.0
    for _ in range(60):
        midpoint = 0.5 * (lower + upper)
        candidate = data_stat + midpoint * mc_cross
        if minimum_covariance_eigenvalue(candidate, component_labels) >= -tolerance:
            lower = midpoint
        else:
            upper = midpoint
    covariance = data_stat + lower * mc_cross
    return lower, 0.5 * (covariance + covariance.T)


# Find exact connected components in one union of covariance supports
def covariance_component_labels(*matrices: np.ndarray, vectors: np.ndarray | None = None) -> np.ndarray:
    if not matrices:
        return np.empty(0, dtype=np.int64)
    size = matrices[0].shape[0]
    graph = sparse.csr_matrix((size, size), dtype=np.int8)
    for matrix in matrices:
        values = sparse.csr_matrix(matrix)
        if values.shape != (size, size):
            raise ValueError("Covariance component matrices have incompatible dimensions")
        support = (values != 0.0).astype(np.int8)
        support.setdiag(0)
        support.eliminate_zeros()
        graph += support
    if vectors is not None:
        for vector in np.asarray(vectors).T:
            indices = np.flatnonzero(vector)
            if len(indices) > 1:
                graph += sparse.csr_matrix((np.ones(len(indices), dtype=np.int8),
                                            (np.full(len(indices), indices[0]), indices)), shape=(size, size))
    graph.data[:] = 1
    _, labels = sparse.csgraph.connected_components(graph, directed=False, return_labels=True)
    return np.asarray(labels, dtype=np.int64)


# Group retained-bin indices by their exact covariance component
def component_indices(labels):
    order = np.argsort(labels, kind="stable")
    return [] if not len(order) else np.split(order, np.flatnonzero(np.diff(np.asarray(labels)[order])) + 1)


# Encode one block covariance and its collective nuisance vectors for persistence
def covariance_decomposition(
    block_covariance: np.ndarray, collective_vectors: np.ndarray, *, component_labels: np.ndarray | None = None
) -> dict:
    block_covariance = sparse.csr_matrix(block_covariance, dtype=float)
    collective_vectors = np.asarray(collective_vectors, dtype=float)
    if block_covariance.ndim != 2 or block_covariance.shape[0] != block_covariance.shape[1]:
        raise ValueError("Block covariance must be square")
    size = block_covariance.shape[0]
    if collective_vectors.ndim != 2 or collective_vectors.shape[0] != size:
        raise ValueError("Collective covariance vectors have incompatible dimensions")
    labels = (
        covariance_component_labels(block_covariance)
        if component_labels is None
        else np.asarray(component_labels, dtype=np.int64)
    )
    if labels.shape != (size,):
        raise ValueError("Covariance component labels have incompatible dimensions")
    matrix = sparse.csr_matrix(block_covariance)
    matrix.eliminate_zeros()
    return {
        "data_block_values": np.asarray(matrix.data, dtype=float),
        "data_block_columns": np.asarray(matrix.indices, dtype=np.int64),
        "data_block_row_offsets": np.asarray(matrix.indptr, dtype=np.int64),
        "data_block_shape": np.asarray(matrix.shape, dtype=np.int64),
        "data_collective_vectors": collective_vectors,
        "data_component_labels": labels,
    }


# Reconstruct and cache dense covariance blocks used by trial likelihoods
def prepare_covariance_structure(payload: dict) -> dict | None:
    cached = payload.get("_covariance_structure")
    if cached is not None:
        return cached
    if any(key not in payload for key in STRUCTURED_COVARIANCE_KEYS):
        return None
    shape = tuple(np.asarray(payload["data_block_shape"], dtype=np.int64).tolist())
    if len(shape) != 2 or shape[0] != shape[1]:
        raise ValueError("Stored block covariance shape is invalid")
    matrix = sparse.csr_matrix(
        (
            np.asarray(payload["data_block_values"], dtype=float),
            np.asarray(payload["data_block_columns"], dtype=np.int64),
            np.asarray(payload["data_block_row_offsets"], dtype=np.int64),
        ),
        shape=shape,
    )
    matrix.sort_indices()
    labels = np.asarray(payload["data_component_labels"], dtype=np.int64)
    vectors = np.asarray(payload["data_collective_vectors"], dtype=float)
    if labels.shape != (shape[0],) or vectors.ndim != 2 or vectors.shape[0] != shape[0]:
        raise ValueError("Stored covariance decomposition dimensions are invalid")
    rows, columns = matrix.nonzero()
    if np.any(labels[rows] != labels[columns]):
        raise ValueError("Stored covariance components do not contain the block support")
    components = []
    for indices in component_indices(labels):
        components.append(
            {"indices": indices, "covariance": np.asarray(matrix[indices][:, indices].toarray(), dtype=float)}
        )
    structure = {
        "size": shape[0],
        "matrix": matrix,
        "components": tuple(components),
        "collective_vectors": vectors,
        "diagonal": np.asarray(matrix.diagonal(), dtype=float),
    }
    payload["_covariance_structure"] = structure
    return structure


# Select covariance coordinates while retaining their independent block structure
def select_covariance_structure(structure: dict, indices: np.ndarray) -> dict:
    indices = np.asarray(indices, dtype=np.int64)
    size = int(structure["size"])
    if indices.ndim != 1 or np.any((indices < 0) | (indices >= size)):
        raise ValueError("Structured covariance selection indices are invalid")
    if len(np.unique(indices)) != len(indices):
        raise ValueError("Structured covariance selection indices are not unique")
    position = np.full(size, -1, dtype=np.int64)
    position[indices] = np.arange(len(indices), dtype=np.int64)
    retained = position >= 0
    components = []
    for component in structure["components"]:
        component_indices = np.asarray(component["indices"], dtype=np.int64)
        selected = retained[component_indices]
        if not np.any(selected):
            continue
        components.append(
            {
                "indices": position[component_indices[selected]],
                "covariance": np.asarray(component["covariance"], dtype=float)[np.ix_(selected, selected)],
            }
        )
    vectors = np.asarray(structure["collective_vectors"], dtype=float)[indices]
    diagonal = np.asarray(structure["diagonal"], dtype=float)[indices]
    return {"size": len(indices), "components": tuple(components), "collective_vectors": vectors, "diagonal": diagonal}


# Build one fixed data covariance payload from data and MC correlation inputs
def build_data_covariance(
    *, data: list, obs: list, mcdata_by_dataset: list, mc_correlation_datacards: list[dict] | None = None
) -> dict:
    layout = build_layout(data)
    terms = _data_covariance_terms(data, layout, sparse_output=True)
    data_stat_local = terms["data_stat_local"]
    data_syst_local = terms["data_syst_local"]
    data_count_stat_uncertainty = terms["data_count_stat_uncertainty"]
    collective_vectors = terms["data_collective_vectors"]
    collective_sources = terms["data_collective_sources"]
    statistical_vectors = collective_vectors[
        :, [source.get("category") == "statistical" for source in collective_sources]
    ]
    data_stat_uncertainty = np.sqrt(np.clip(data_stat_local.diagonal() + np.sum(statistical_vectors**2, axis=1), 0.0, None))
    mc_cross_correlation = estimate_mc_cross_correlation(
        obs=obs, mcdata_by_dataset=mcdata_by_dataset, layout=layout, sparse_output=True)
    scaling = sparse.diags(data_count_stat_uncertainty)
    mc_cross_covariance = scaling @ mc_cross_correlation @ scaling
    component_labels = covariance_component_labels(data_stat_local, data_syst_local, mc_cross_covariance)
    mc_correlation_scale = sparse_mc_covariance_scale(data_stat_local, mc_cross_covariance, statistical_vectors)
    data_block_covariance = data_stat_local + data_syst_local + mc_correlation_scale * mc_cross_covariance
    payload = CovariancePayload({
        "schema_version": DATA_COVARIANCE_SCHEMA,
        "method": "block_data_covariance_with_collective_sources",
        "asymmetric_rule": "mean_absolute_up_down",
        "layout": layout,
        "data_uncertainty_sources": terms["data_uncertainty_sources"],
        "data_collective_sources": collective_sources,
        "mc_correlation_samples": mc_correlation_sample_metadata(mcdata_by_dataset, mc_correlation_datacards),
        "mc_correlation_scale": mc_correlation_scale,
        "data_stat_uncertainty": data_stat_uncertainty,
        "data_count_stat_uncertainty": data_count_stat_uncertainty,
        **sparse_covariance_arrays(data_stat_local, "stat_local"),
        **sparse_covariance_arrays(data_syst_local, "syst_local"),
        **sparse_covariance_arrays(mc_cross_correlation, "mc_cross"),
        **covariance_decomposition(data_block_covariance, collective_vectors, component_labels=component_labels),
    })
    prepare_covariance_structure(payload)
    return payload


# Apply the same finite-sample PSD correction over independent statistical components
def sparse_mc_covariance_scale(statistical, cross, vectors):
    labels = covariance_component_labels(statistical, cross, vectors=vectors)
    factor = 1.0
    for indices in component_indices(labels):
        block = statistical[indices][:, indices].toarray()
        selected = vectors[indices]
        block += selected @ selected.T
        local_factor, _ = scale_mc_cross_covariance_to_psd(
            block, cross[indices][:, indices].toarray())
        factor = min(factor, local_factor)
    return factor


# Save the fixed data covariance matrices and readable metadata
def save_data_covariance(payload: dict, *, output_dir: str) -> tuple[str, str]:
    validate_data_covariance(payload)
    ensure_dir(output_dir, exist_ok=True)
    npz_path = os.path.join(output_dir, "data_covariance.npz")
    json_path = os.path.join(output_dir, "data_covariance.json")
    npz_tmp = f"{npz_path}.tmp"
    json_tmp = f"{json_path}.tmp"
    with open(npz_tmp, "wb") as output:
        np.savez_compressed(output, **{key: payload[key] for key in DATA_COVARIANCE_ARRAYS})
        output.flush()
        os.fsync(output.fileno())
    metadata = {
        key: value for key, value in payload.items() if key not in DATA_COVARIANCE_ARRAYS and not key.startswith("_")
    }
    metadata["matrix_file"] = os.path.basename(npz_path)
    metadata["arrays"] = list(DATA_COVARIANCE_ARRAYS)
    metadata["matrix_sha256"] = file_sha256(npz_tmp)
    with open(json_tmp, "w", encoding="utf-8") as output:
        json.dump(metadata, output, indent=2, sort_keys=True)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    os.replace(npz_tmp, npz_path)
    os.replace(json_tmp, json_path)
    return npz_path, json_path


# Load and integrity-check one fixed data covariance payload
def load_data_covariance(*, npz_path: str, json_path: str, data=None) -> dict:
    payload = CovariancePayload(load_json_file(json_path))
    if payload.get("matrix_file") != os.path.basename(npz_path):
        raise ValueError("Data covariance matrix filename does not match metadata")
    if payload.get("matrix_sha256") != file_sha256(npz_path):
        raise ValueError("Data covariance matrix checksum does not match metadata")
    if payload.get("arrays") != list(DATA_COVARIANCE_ARRAYS):
        raise ValueError("Data covariance array manifest is invalid")
    with np.load(npz_path) as arrays:
        if sorted(arrays.files) != sorted(DATA_COVARIANCE_ARRAYS):
            raise ValueError("Data covariance NPZ contents do not match metadata")
        for key in DATA_COVARIANCE_ARRAYS:
            payload[key] = np.asarray(arrays[key])
    validate_data_covariance(payload, data=data)
    prepare_covariance_structure(payload)
    return payload


# Check a dense covariance against its sparse block and low-rank decomposition
def validate_covariance_decomposition(covariance: np.ndarray, structure: dict) -> None:
    covariance = np.asarray(covariance, dtype=float)
    size = int(structure["size"])
    if covariance.shape != (size, size):
        raise ValueError("Structured covariance has incompatible dimensions")
    matrix = structure["matrix"]
    vectors = structure["collective_vectors"]
    for start in range(0, size, 256):
        stop = min(start + 256, size)
        expected = np.asarray(matrix[start:stop].toarray(), dtype=float)
        expected += vectors[start:stop] @ vectors.T
        if not np.allclose(covariance[start:stop], expected, rtol=1.0e-10, atol=1.0e-14):
            raise ValueError("Data covariance structured decomposition is inconsistent")


# Validate one covariance after removing its collective and optional fixed terms
def validate_local_covariance(
    covariance: np.ndarray,
    collective_vectors: np.ndarray,
    component_labels: np.ndarray,
    *,
    name: str,
    subtract: np.ndarray | None = None,
) -> None:
    covariance = np.asarray(covariance, dtype=float)
    vectors = np.asarray(collective_vectors, dtype=float)
    labels = np.asarray(component_labels, dtype=np.int64)
    removed = None if subtract is None else np.asarray(subtract, dtype=float)
    size = len(labels)
    for start in range(0, size, 256):
        stop = min(start + 256, size)
        collective = vectors[start:stop] @ vectors.T
        local = covariance[start:stop] - collective
        if removed is not None:
            local -= removed[start:stop]
        outside = labels[start:stop, None] != labels[None, :]
        reference = np.abs(covariance[start:stop]) + np.abs(collective)
        if removed is not None:
            reference += np.abs(removed[start:stop])
        tolerance = 1.0e-14 + 1.0e-10 * reference
        if np.any(np.abs(local[outside]) > tolerance[outside]):
            raise ValueError(f"{name} crosses independent covariance components")
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        selector = np.ix_(indices, indices)
        block_vectors = vectors[indices]
        collective = block_vectors @ block_vectors.T
        block = covariance[selector] - collective
        if removed is not None:
            block -= removed[selector]
        reference = max(float(np.max(np.abs(covariance[selector]))), float(np.max(np.abs(collective))), 1.0e-30)
        if removed is not None:
            reference = max(reference, float(np.max(np.abs(removed[selector]))))
        if float(np.min(np.linalg.eigvalsh(block))) < -1.0e-10 * reference:
            raise ValueError(f"{name} is not positive semidefinite")


# Validate sparse covariance terms and their independent physical components
def validate_data_covariance(payload: dict, *, data: list | None = None) -> None:
    if payload.get("schema_version") != DATA_COVARIANCE_SCHEMA:
        raise ValueError("Unsupported data covariance schema")
    layout = payload.get("layout")
    if not isinstance(layout, list):
        raise ValueError("Data covariance layout must be a list")
    size = 0 if not layout else int(layout[-1]["stop"])
    structure = prepare_covariance_structure(payload)
    if structure is None or structure["size"] != size:
        raise ValueError("Data covariance decomposition does not match the layout")
    uncertainty = np.asarray(payload["data_stat_uncertainty"], dtype=float)
    count_uncertainty = np.asarray(payload["data_count_stat_uncertainty"], dtype=float)
    vectors = np.asarray(structure["collective_vectors"], dtype=float)
    sources = payload.get("data_collective_sources")
    labels = np.asarray(payload["data_component_labels"], dtype=np.int64)
    if not isinstance(sources, list) or vectors.shape != (size, len(sources)):
        raise ValueError("Data collective covariance source dimensions are invalid")
    if labels.shape != (size,) or np.any(labels < 0):
        raise ValueError("Data covariance component labels are invalid")
    if not np.all(np.isfinite(vectors)):
        raise ValueError("Data collective covariance vectors contain non-finite values")
    for values in (uncertainty, count_uncertainty):
        if values.shape != (size,) or not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("Data statistical uncertainties are invalid")
    if np.any(count_uncertainty > uncertainty + 1.0e-12 * np.maximum(uncertainty, 1.0)):
        raise ValueError("Data counting uncertainties exceed total statistical uncertainties")
    matrices = {name: sparse_covariance_matrix(payload, name) for name in ("stat_local", "syst_local", "mc_cross")}
    for name, matrix in matrices.items():
        matrix.check_format(full_check=True)
        if matrix.shape != (size, size) or not np.all(np.isfinite(matrix.data)):
            raise ValueError(f"Invalid sparse covariance: {name}")
        difference = (matrix - matrix.T).tocoo()
        if difference.nnz and not np.allclose(difference.data, 0.0, rtol=1.0e-10, atol=1.0e-14):
            raise ValueError("Data covariance matrices are not symmetric")
        rows, columns = matrix.nonzero()
        if np.any(labels[rows] != labels[columns]):
            raise ValueError(f"{name} crosses independent covariance components")
    cross = matrices["mc_cross"]
    if np.any(np.abs(cross.data) > 1.0 + 1.0e-10):
        raise ValueError("MC cross-correlation is outside [-1, 1]")
    offset = 0
    for record in layout:
        if record["start"] != offset or record["stop"] - offset != len(record["bins"]):
            raise ValueError("Data covariance layout has inconsistent offsets")
        indices = record_indices(record)
        block = cross[indices][:, indices]
        if block.nnz and not np.allclose(block.data, 0.0, rtol=0.0, atol=1.0e-14):
            raise ValueError("MC cross-correlation has same-observable entries")
        offset = record["stop"]
    factor = float(payload["mc_correlation_scale"])
    if not np.isfinite(factor) or not 0.0 <= factor <= 1.0:
        raise ValueError("MC correlation scale is outside [0, 1]")
    scaling = sparse.diags(count_uncertainty)
    corrected = matrices["stat_local"] + factor * (scaling @ cross @ scaling)
    expected = corrected + matrices["syst_local"]
    difference = (expected - structure["matrix"]).tocoo()
    if difference.nnz and not np.allclose(difference.data, 0.0, rtol=1.0e-10, atol=1.0e-14):
        raise ValueError("Data covariance structured decomposition is inconsistent")
    statistical = np.asarray([source["category"] == "statistical" for source in sources], dtype=bool)
    diagonal = corrected.diagonal() + np.sum(vectors[:, statistical] ** 2, axis=1)
    if not np.allclose(diagonal, uncertainty**2, rtol=1.0e-10, atol=1.0e-14):
        raise ValueError("Data statistical covariance diagonal is inconsistent")
    for indices in component_indices(labels):
        for matrix in (matrices["stat_local"], matrices["syst_local"], corrected):
            block = matrix[indices][:, indices].toarray()
            reference = max(float(np.max(np.abs(block), initial=0.0)), 1.0e-30)
            if minimum_covariance_eigenvalue(block) < -1.0e-10 * reference:
                raise ValueError("Data covariance component is not positive semidefinite")
    if data is not None and layout != build_layout(data):
        raise ValueError("Data covariance layout does not match initialized data")


# Flatten trial MC and data vectors using the fixed covariance layout
def comparison_vectors(
    results: dict, layout: list[dict]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vectors = {
        name: [] for name in ("mc_prediction", "data_values", "mc_stat_uncertainty", "data_uncertainty", "fit_weight")
    }
    for record in layout:
        bins = np.asarray(record["bins"], dtype=np.int64)
        mc_item = results["mc"][record["dataset"]][record["subset"]][record["observable"]]
        data_item = results["data"][record["dataset"]][record["subset"]][record["observable"]]
        values = (
            mc_item["hdata"].counts_scaled[bins],
            np.asarray(data_item["hdata"].counts_scaled)[bins],
            mc_item["hdata"].errs_scaled[bins],
            np.asarray(data_item["hdata"].errs_scaled)[bins],
            np.full(len(bins), float(data_item["fitw"])),
        )
        for name, value in zip(vectors, values, strict=True):
            vectors[name].append(value)
    return tuple(array.concatenate(values) for values in vectors.values())


# Flatten the trial MC validity mask using the fixed data layout
def comparison_mc_valid(results: dict, layout: list[dict]) -> np.ndarray:
    valid = []
    for record in layout:
        histogram = results["mc"][record["dataset"]][record["subset"]][record["observable"]]["hdata"]
        valid.extend(np.asarray(histogram.valid, dtype=bool)[record["bins"]])
    return np.asarray(valid, dtype=bool)


# Map full flattened histogram arrays into fixed covariance ordering
def manifest_layout_indices(layout: list[dict], manifest: list[dict]) -> np.ndarray:
    records = {(int(item["dataset"]), int(item["subset"]), str(item["observable"])): item for item in manifest}
    indices = []
    for record in layout:
        key = (int(record["dataset"]), int(record["subset"]), str(record["observable"]))
        if key not in records:
            raise ValueError(f"Full histogram manifest is missing covariance record {key}")
        item = records[key]
        width = int(item["stop"]) - int(item["start"])
        bins = np.asarray(record["bins"], dtype=np.int64)
        if np.any((bins < 0) | (bins >= width)):
            raise ValueError(f"Data covariance bins exceed histogram record {key}")
        indices.extend(int(item["start"]) + bins)
    return np.asarray(indices, dtype=np.int64)


# Expand per-observable fit weights into full flattened bin ordering
def manifest_fit_weights(manifest: list[dict]) -> np.ndarray:
    if not manifest:
        return np.empty(0, dtype=float)
    size = int(manifest[-1]["stop"])
    weights = np.zeros(size, dtype=float)
    for item in manifest:
        weights[int(item["start"]) : int(item["stop"])] = float(item.get("fitw", 1.0))
    return weights


# Compute whether one simulation payload requests the full data covariance
def uses_full_data_covariance(payload: dict) -> bool:
    return payload.get("simdriver") == "GRANIITTI" and payload.get("data_covariance_mode") == "full"


# Load and align one run covariance to positive-weight manifest bins
def load_manifest_covariance(
    *, cdir: str, run_name: str, manifest: list[dict], fit_weights: np.ndarray
) -> tuple[dict, np.ndarray, np.ndarray]:
    covariance_root = os.path.join(cdir, "runs", "icetune", run_name, "results")
    payload = load_data_covariance(
        npz_path=os.path.join(covariance_root, "data_covariance.npz"),
        json_path=os.path.join(covariance_root, "data_covariance.json"),
    )
    indices = manifest_layout_indices(payload["layout"], manifest)
    retained = np.isfinite(fit_weights[indices]) & (fit_weights[indices] > 0.0)
    if not np.any(retained):
        raise ValueError("Data covariance has no positive-weight histogram bins")
    covariance = covariance_selection(payload, indices=np.flatnonzero(retained))
    structure = prepare_covariance_structure(payload)
    if structure is not None:
        payload["_selected_covariance_structure"] = select_covariance_structure(structure, np.flatnonzero(retained))
    return payload, indices[retained], covariance


# Compute one marginal data covariance block in local histogram ordering
def observable_covariance_block(
    *, payload: dict, dataset: int, subset: int, observable: str, bin_count: int
) -> np.ndarray:
    bins, covariance = observable_covariance_selection(
        payload=payload, dataset=dataset, subset=subset, observable=observable, bin_count=bin_count
    )
    block = np.zeros((bin_count, bin_count), dtype=float)
    block[np.ix_(bins, bins)] = covariance
    return block


# Select one observable and its marginal covariance in retained-bin ordering
def observable_covariance_selection(
    *, payload: dict, dataset: int, subset: int, observable: str, bin_count: int
) -> tuple[np.ndarray, np.ndarray]:
    matches = [
        record
        for record in payload["layout"]
        if int(record["dataset"]) == int(dataset)
        and int(record["subset"]) == int(subset)
        and str(record["observable"]) == str(observable)
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one covariance layout record for {(dataset, subset, observable)}")
    record = matches[0]
    bins = np.asarray(record["bins"], dtype=np.int64)
    if np.any((bins < 0) | (bins >= int(bin_count))):
        raise ValueError(f"Data covariance bins exceed histogram {(dataset, subset, observable)}")
    indices = record_indices(record)
    return bins, covariance_selection(payload, indices=indices)


# Assemble MC covariance in the same fixed bin ordering as the data
def comparison_covariance(results, layout):
    histograms = [results["mc"][r["dataset"]][r["subset"]][r["observable"]]["hdata"] for r in layout]
    if all(not h.mc_events and not (h.density and h.density_uncertainty == "shape") for h in histograms):
        return None
    return uncertainty.joint_covariance(histograms, [np.asarray(r["bins"], dtype=int) for r in layout])
