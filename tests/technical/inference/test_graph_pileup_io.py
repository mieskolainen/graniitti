# Tests for graph based pileup baseline data and reconstruction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import pytest
import uproot

RUN_STUDY_DIR = (
    Path(__file__).resolve().parents[3]
    / "tests"
    / "physics"
    / "studies"
    / "pileup_kron"
)
sys.path.insert(0, str(RUN_STUDY_DIR))

import analyze_delphes as baseline  # noqa: E402


# Compute a recovery argument namespace matching the analyzer defaults
def make_recovery_args(**overrides: float | int) -> argparse.Namespace:
    values: dict[str, float | int] = {
        "graph_kron_recovery_min_hard_fraction": 0.20,
        "graph_kron_recovery_min_hard_pt": 15.0,
        "graph_kron_recovery_pt_min": 15.0,
        "graph_kron_recovery_dedupe_radius": 0.28,
        "graph_kron_recovery_max_add": 0,
        "graph_kron_recovery_residual_fraction": 0.90,
        "graph_kron_recovery_score_pileup_penalty": 0.50,
        "graph_kron_recovery_log_odds_min": 2.0,
        "graph_kron_recovery_likelihood_weight": 1.0,
        "graph_kron_recovery_output_scale": 1.0,
        "graph_kron_latent_recovery_raw_pt_min": 70.0,
        "graph_kron_latent_recovery_min_hard_pt": 1.0,
        "graph_kron_latent_recovery_max_add": 1,
        "graph_kron_latent_recovery_output_scale": 0.30,
        "graph_kron_neutral_pileup_scale": 1.0,
        "graph_kron_subtracted_hard_strength": 0.0,
        "reco_pt_min": 15.0,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


# Decode real ROOT baskets across precision, empty events and basket boundaries
@pytest.mark.parametrize('dtype', ['float32', 'float64', 'int32'])
def test_numeric_branch_round_trip(tmp_path, dtype):
    path = tmp_path / 'delphes.root'
    expected = ak.to_list(ak.values_astype(ak.Array([[1.25, 2.5], [], [3.75], [4.5, 5.25, 6.5 + 2**-30]]), dtype))
    with uproot.recreate(path) as output:
        tree = output.mktree('Delphes', {'Jet.PT': f'var * {dtype}'}, counter_name=lambda _: 'Jet')
        for events in (expected[:2], expected[2:]):
            tree.extend({'Jet.PT': ak.values_astype(ak.Array(events), dtype)})
    with uproot.open(path) as source:
        tree = source['Delphes']
        decode = baseline.load_jagged_int_branch if dtype == 'int32' else baseline.load_jagged_float_branch
        assert decode(tree, 'Jet.PT') == expected
        assert baseline.load_count_branch(tree, 'Jet') == [2, 0, 1, 3]
        assert tree['Jet.PT'].num_baskets == 2


# Scalar ROOT baskets have no jagged offsets and exercise the count fallback directly
def test_offsets_can_fall_back_to_counts(tmp_path):
    path = tmp_path / 'counts.root'
    with uproot.recreate(path) as output:
        tree = output.mktree('Delphes', {'Jet': 'int32'})
        tree.extend({'Jet': np.array([1, 2], dtype=np.int32)})
    with uproot.open(path) as source:
        basket = source['Delphes']['Jet'].basket(0)
        assert basket.byte_offsets is None
        offsets = baseline.basket_item_offsets(basket, 4, [1, 2], 0, 2)
        assert offsets == [0, 1, 3]
        baseline.validate_item_offsets('Jet.PT', offsets, 3, 0, 2)
        with pytest.raises(ValueError, match='expected 4'):
            baseline.validate_item_offsets('Jet.PT', offsets, 3, 0, 3)


# Check that branch length mismatches fail before physics metrics are computed
def test_jagged_length_mismatch() -> None:
    with pytest.raises(ValueError, match="branch length mismatch"):
        baseline.validate_equal_event_lengths(
            "Jet",
            {
                "Jet.PT": [[1.0, 2.0]],
                "Jet.Eta": [[0.1]],
            },
        )


# Check that invalid count branches are rejected
def test_negative_counts_are_rejected(tmp_path):
    path = tmp_path / 'negative.root'
    with uproot.recreate(path) as output:
        tree = output.mktree('Delphes', {'Jet': 'int32'})
        tree.extend({'Jet': np.array([1, -1], dtype=np.int32)})
    with uproot.open(path) as source, pytest.raises(ValueError, match='negative'):
        baseline.load_count_branch(source['Delphes'], 'Jet')


# Check that event-count validation does not require equal per-event object counts
def test_event_collection_sizes() -> None:
    result = baseline.validate_equal_event_counts(
        "Collections",
        {
            "GenJet": [[object()], [object(), object()]],
            "Jet": [[object(), object()], []],
        },
    )
    assert result is None


# Check that zero-reference events do not dilute match-efficiency means
def test_efficiency_zero_reference() -> None:
    matches = [
        baseline.Match(
            event=1,
            collection="Reco",
            gen_pt=50.0,
            gen_eta=0.0,
            gen_phi=0.0,
            reco_pt=48.0,
            reco_eta=0.01,
            reco_phi=0.02,
            dr=0.03,
        )
    ]
    row = baseline.summarize_collection(
        "Reco", matches, gen_counts=[0, 2], reco_counts=[1, 1], fake_counts=[1, 0]
    )
    assert row["match_efficiency_mean"] == pytest.approx(0.5)
    assert row["match_efficiency_global"] == pytest.approx(0.5)


# Check that event slicing keeps the requested deterministic subset
def test_limit_events_uses_leading_slice() -> None:
    events = [[0], [1], [2]]
    assert baseline.limit_events(events, 2) == [[0], [1]]
    assert baseline.limit_events(events, 0) == events
    assert baseline.limit_events(events, 2, 1) == [[1], [2]]
    assert baseline.limit_events(events, 0, 0, [2, 0]) == [[2], [0]]


# Check that event-index parsing accepts inline and file inputs
def test_parse_event_indices(tmp_path: Path) -> None:
    assert baseline.parse_event_indices(
        argparse.Namespace(event_indices="2, 5 7", event_indices_file="")
    ) == [2, 5, 7]
    index_file = tmp_path / "events.txt"
    index_file.write_text("3\n4,8\n", encoding="utf-8")
    assert baseline.parse_event_indices(
        argparse.Namespace(event_indices="", event_indices_file=str(index_file))
    ) == [3, 4, 8]
    with pytest.raises(ValueError, match="unique"):
        baseline.parse_event_indices(argparse.Namespace(event_indices="1,1", event_indices_file=""))


# Check that tiled anti-kT merges particles inside the jet radius
def test_tiled_antikt_merges_close_particles() -> None:
    particles = [
        baseline.WeightedParticle(pt=10.0, eta=0.0, phi=0.0, mass=0.0),
        baseline.WeightedParticle(pt=5.0, eta=0.0, phi=0.0, mass=0.0),
    ]
    jets = baseline.cluster_jets(particles, radius=0.4, pt_min=0.0)
    assert len(jets) == 1
    assert jets[0].pt == pytest.approx(15.0)


# Check that tiled anti-kT keeps separated particles as different jets
def test_tiled_antikt_separated_particles() -> None:
    particles = [
        baseline.WeightedParticle(pt=10.0, eta=0.0, phi=0.0, mass=0.0),
        baseline.WeightedParticle(pt=5.0, eta=1.0, phi=0.0, mass=0.0),
    ]
    jets = baseline.cluster_jets(particles, radius=0.4, pt_min=0.0)
    assert [jet.pt for jet in jets] == [pytest.approx(10.0), pytest.approx(5.0)]


# Check that graph-label propagation preserves simplex normalization
def test_label_propagation_rows_sum_to_one() -> None:
    adjacency = [[(1, 1.0)], [(0, 1.0)]]
    unary_support = np.asarray([[1.0, 0.0], [0.0, 1.0]])
    strengths = np.asarray([5.0, 5.0])
    labels = baseline.solve_label_propagation(adjacency, unary_support, strengths, iterations=8)
    assert np.allclose(np.sum(labels, axis=1), 1.0)
    assert labels[0, 0] > labels[0, 1]
    assert labels[1, 1] > labels[1, 0]


# Check that an isolated neutral with only sink support stays assigned to sink
def test_isolated_neutral_stays_sink_terminal() -> None:
    unary_support = np.asarray([[0.0, 1.0]])
    strengths = np.asarray([1.0])
    labels = baseline.solve_label_propagation([[]], unary_support, strengths, iterations=4)
    assert labels[0, 1] == pytest.approx(1.0)


# Check that neutral local support favours the nearby primary charged anchor
def test_terminal_nearby_primary() -> None:
    candidates = [
        baseline.Candidate(
            pt=10.0, eta=0.0, phi=0.0, mass=0.0, charge=1, vertex=10, kind="track", index=0
        ),
        baseline.Candidate(
            pt=10.0, eta=2.0, phi=0.0, mass=0.0, charge=1, vertex=20, kind="track", index=1
        ),
        baseline.Candidate(
            pt=5.0, eta=0.05, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=2
        ),
    ]
    unary_support, strengths = baseline.build_terminal_support(
        candidates,
        terminals=[10, 20],
        radius=0.4,
        sigma=0.2,
        anchor_strength=20.0,
        local_support_strength=1.0,
        sink_strength=0.1,
    )
    assert unary_support[2, 0] > unary_support[2, 1]
    assert unary_support[2, 0] > unary_support[2, 2]
    assert strengths[0] == pytest.approx(20.0)


# Check that fitted Fisher odds separate charged hard and pile-up patterns
def test_fisher_model_separates_charged_patterns() -> None:
    hard_features = [
        baseline.fisher_feature_vector(
            baseline.Candidate(
                pt=12.0 + index,
                eta=0.1,
                phi=0.0,
                mass=0.0,
                charge=1,
                vertex=0,
                kind="track",
                index=index,
            ),
            primary_support=8.0 + index,
            pileup_support=0.2,
        )
        for index in range(6)
    ]
    pileup_features = [
        baseline.fisher_feature_vector(
            baseline.Candidate(
                pt=2.0 + index,
                eta=1.5,
                phi=0.0,
                mass=0.0,
                charge=1,
                vertex=1,
                kind="track",
                index=index,
            ),
            primary_support=0.2,
            pileup_support=8.0 + index,
        )
        for index in range(6)
    ]
    model = baseline.fit_fisher_model(
        hard_features, pileup_features, shrinkage=0.5, min_class_count=2
    )
    assert model is not None
    assert baseline.fisher_log_odds(model, hard_features[0], model.sample_log_prior_odds) > 0.0
    assert baseline.fisher_log_odds(model, pileup_features[0], model.sample_log_prior_odds) < 0.0


# Check that Fisher likelihood support strengthens a primary-like neutral
def test_terminal_support_fisher_likelihood_neutral() -> None:
    candidates = [
        baseline.Candidate(
            pt=10.0, eta=0.0, phi=0.0, mass=0.0, charge=1, vertex=0, kind="track", index=0
        ),
        baseline.Candidate(
            pt=10.0, eta=1.5, phi=0.0, mass=0.0, charge=1, vertex=1, kind="track", index=1
        ),
        baseline.Candidate(
            pt=5.0, eta=0.05, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=2
        ),
    ]
    model = baseline.FisherGraphModel(
        weights=np.asarray([0.0, 2.0, -2.0, 1.0, 2.0, -2.0, 1.0, 0.0]),
        bias=0.0,
        sample_log_prior_odds=0.0,
        hard_mean=np.zeros(8),
        pileup_mean=np.zeros(8),
        covariance=np.eye(8),
        hard_count=10,
        pileup_count=10,
    )
    local_support, _strengths = baseline.build_terminal_support(
        candidates,
        terminals=[0, 1],
        radius=0.4,
        sigma=0.2,
        anchor_strength=20.0,
        local_support_strength=1.0,
        sink_strength=0.1,
    )
    fisher_support, fisher_strengths = baseline.build_terminal_support(
        candidates,
        terminals=[0, 1],
        radius=0.4,
        sigma=0.2,
        anchor_strength=20.0,
        local_support_strength=1.0,
        sink_strength=0.1,
        fisher_model=model,
        fisher_strength=3.0,
        fisher_event_prior_weight=0.0,
    )
    assert fisher_support[2, 0] > local_support[2, 0]
    assert fisher_strengths[2] == pytest.approx(4.0)


# Check that neutral-orphan support can raise unsupported neutral hard probability
def test_terminal_support_orphan_neutral_prior() -> None:
    candidates = [
        baseline.Candidate(
            pt=20.0, eta=0.0, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=0
        ),
    ]
    baseline_support, _baseline_strengths = baseline.build_terminal_support(
        candidates,
        [0],
        radius=0.4,
        sigma=0.2,
        anchor_strength=10.0,
        local_support_strength=1.0,
        sink_strength=0.5,
    )
    orphan_support, _orphan_strengths = baseline.build_terminal_support(
        candidates,
        [0],
        radius=0.4,
        sigma=0.2,
        anchor_strength=10.0,
        local_support_strength=1.0,
        sink_strength=0.5,
        orphan_neutral_strength=3.0,
        orphan_neutral_pt_scale=5.0,
        orphan_neutral_pileup_scale=1.0,
    )
    assert orphan_support[0, 0] > baseline_support[0, 0]


# Check that long-range radiation edges connect charged labels to neutrals
def test_radiation_edges_connect_charged_neutral() -> None:
    candidates = [
        baseline.Candidate(
            pt=20.0, eta=0.0, phi=0.0, mass=0.0, charge=1, vertex=0, kind="track", index=0
        ),
        baseline.Candidate(
            pt=5.0, eta=0.8, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=1
        ),
    ]
    local_only = baseline.build_candidate_graph(
        candidates,
        radius=0.4,
        sigma=0.2,
        radiation_radius=1.0,
        radiation_core=0.04,
        radiation_angular_power=1.0,
        radiation_pt_power=1.0,
        radiation_edge_strength=0.0,
        radiation_edge_max_neighbours=16,
    )
    with_radiation = baseline.build_candidate_graph(
        candidates,
        radius=0.4,
        sigma=0.2,
        radiation_radius=1.0,
        radiation_core=0.04,
        radiation_angular_power=1.0,
        radiation_pt_power=1.0,
        radiation_edge_strength=0.2,
        radiation_edge_max_neighbours=16,
    )
    assert local_only[1] == []
    assert any(index == 0 and weight > 0.0 for index, weight in with_radiation[1])


# Check that the compiled CSR graph solver preserves label normalization
def test_numba_csr_graph_solver_labels() -> None:
    candidates = [
        baseline.Candidate(
            pt=20.0, eta=0.0, phi=0.0, mass=0.0, charge=1, vertex=0, kind="track", index=0
        ),
        baseline.Candidate(
            pt=5.0, eta=0.2, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=1
        ),
        baseline.Candidate(
            pt=10.0, eta=1.0, phi=0.0, mass=0.0, charge=1, vertex=1, kind="track", index=2
        ),
    ]
    indptr, indices, weights = baseline.build_candidate_graph_csr(
        candidates,
        radius=0.4,
        sigma=0.2,
        radiation_radius=1.2,
        radiation_core=0.04,
        radiation_angular_power=1.0,
        radiation_pt_power=1.0,
        radiation_edge_strength=0.1,
        radiation_edge_max_neighbours=4,
    )
    unary_support = np.asarray([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
    strengths = np.asarray([20.0, 1.0, 20.0])
    labels = baseline.solve_label_propagation_csr_numba(
        indptr, indices, weights, unary_support, strengths, 4
    )
    assert np.allclose(np.sum(labels, axis=1), 1.0)
    assert labels[1, 0] > 0.5


# Check that tracks without Delphes VertexIndex are assigned by nearest vertex z
def test_track_vertex_assignment_nearest_z() -> None:
    vertex = baseline.assign_track_vertex(
        track_z=2.2,
        vertex_indices=[0, 7, 12],
        vertex_z=[-10.0, 2.0, 25.0],
        fallback=-1,
    )
    assert vertex == 7


# Check that graph jet filtering keeps only jets with hard-scatter support
def test_filter_rejects_pileup_supported_jets() -> None:
    jets = [
        baseline.Jet(pt=40.0, eta=0.0, phi=0.0, mass=5.0, index=0),
        baseline.Jet(pt=35.0, eta=2.0, phi=0.0, mass=4.0, index=1),
    ]
    candidates = [
        baseline.Candidate(
            pt=25.0, eta=0.03, phi=0.0, mass=0.0, charge=1, vertex=0, kind="track", index=0
        ),
        baseline.Candidate(
            pt=20.0, eta=2.02, phi=0.0, mass=0.0, charge=1, vertex=3, kind="track", index=1
        ),
    ]
    kept = baseline.filter_graph_supported_jets(
        jets,
        candidates,
        hard_weights=[1.0, 0.0],
        radius=0.4,
        min_fraction=0.5,
        min_hard_pt=10.0,
        residual_fraction=0.0,
    )
    assert len(kept) == 1
    assert kept[0].pt == pytest.approx(40.0)


# Check that weighted graph jets recover residual local candidate energy
def test_graph_residual_calibration() -> None:
    jets = [baseline.Jet(pt=20.0, eta=0.0, phi=0.0, mass=2.0, index=0)]
    candidates = [
        baseline.Candidate(
            pt=20.0, eta=0.01, phi=0.0, mass=0.0, charge=1, vertex=0, kind="track", index=0
        ),
        baseline.Candidate(
            pt=10.0, eta=0.02, phi=0.0, mass=0.0, charge=0, vertex=-1, kind="neutral", index=1
        ),
    ]
    calibrated = baseline.calibrate_weighted_graph_jets(
        jets,
        candidates,
        hard_weights=[1.0, 0.0],
        radius=0.4,
        min_fraction=0.1,
        min_hard_pt=1.0,
        residual_fraction=0.5,
    )
    assert len(calibrated) == 1
    assert calibrated[0].pt == pytest.approx(25.0)
    assert calibrated[0].mass == pytest.approx(2.5)


# Check that candidate limiting preserves neutral energy candidates
def test_candidate_limit_is_not_charged_only() -> None:
    candidates = [
        baseline.Candidate(
            pt=100.0 - index,
            eta=0.0,
            phi=0.0,
            mass=0.0,
            charge=1,
            vertex=index + 1,
            kind="track",
            index=index,
        )
        for index in range(20)
    ]
    candidates.extend(
        baseline.Candidate(
            pt=50.0 - index,
            eta=0.0,
            phi=0.1,
            mass=0.0,
            charge=0,
            vertex=-1,
            kind="neutral",
            index=20 + index,
        )
        for index in range(20)
    )
    limited = baseline.limit_candidates(candidates, max_candidates=10, primary_vertex=0)
    assert len(limited) == 10
    assert any(candidate.charge == 0 for candidate in limited)
    assert any(candidate.charge != 0 for candidate in limited)
    assert [candidate.index for candidate in limited] == list(range(10))


# Check that graph recovery cannot create a jet without raw graph anti-kT support
def test_recovery_requires_raw_graph_candidate() -> None:
    args = make_recovery_args()
    recovered = baseline.append_graph_recovery_jets([], [], args)
    assert recovered == []


# Check that graph recovery uses only graph support for pT
def test_recovery_graph_support_for_pt() -> None:
    args = make_recovery_args(graph_kron_recovery_max_add=3)
    raw_records = [
        baseline.GraphJetRecord(
            jet=baseline.Jet(pt=30.0, eta=0.0, phi=0.0, mass=3.0, index=0),
            hard_fraction=2.0 / 3.0,
            hard_pt=20.0,
            total_pt=30.0,
            log_odds=4.0,
        )
    ]
    recovered = baseline.append_graph_recovery_jets([], raw_records, args)
    expected_pt = 20.0 + 0.90 * 10.0
    assert len(recovered) == 1
    assert recovered[0].pt == pytest.approx(expected_pt)
    assert recovered[0].eta == pytest.approx(0.0)
    assert recovered[0].phi == pytest.approx(0.0)


# Check that charged-subtracted support can recover neutral orphan pT
def test_hard_pt_neutral_subtraction() -> None:
    record = baseline.GraphJetRecord(
        jet=baseline.Jet(pt=50.0, eta=0.0, phi=0.0, mass=5.0, index=0),
        hard_fraction=0.1,
        hard_pt=5.0,
        total_pt=50.0,
        primary_charged_pt=4.0,
        pileup_charged_pt=6.0,
        neutral_pt=40.0,
        neutral_hard_pt=1.0,
    )
    assert baseline.record_effective_hard_pt(
        record, neutral_pileup_scale=1.0, subtracted_strength=0.0
    ) == pytest.approx(5.0)
    assert baseline.record_effective_hard_pt(
        record, neutral_pileup_scale=1.0, subtracted_strength=1.0
    ) == pytest.approx(38.0)


# Check that graph recovery does not duplicate an existing graph base jet
def test_recovery_deduplicates_base_jets() -> None:
    args = make_recovery_args(graph_kron_recovery_max_add=3)
    base = [baseline.Jet(pt=80.0, eta=0.1, phi=0.2, mass=8.0, index=0)]
    raw_records = [
        baseline.GraphJetRecord(
            jet=baseline.Jet(pt=30.0, eta=0.1, phi=0.2, mass=3.0, index=0),
            hard_fraction=1.0,
            hard_pt=30.0,
            total_pt=30.0,
            log_odds=4.0,
        )
    ]
    recovered = baseline.append_graph_recovery_jets(base, raw_records, args)
    assert len(recovered) == 1
    assert recovered[0].pt == pytest.approx(80.0)


# Check that recovery rejects candidates with low graph likelihood odds
def test_recovery_rejects_low_likelihood_odds() -> None:
    args = make_recovery_args(graph_kron_recovery_max_add=3)
    raw_records = [
        baseline.GraphJetRecord(
            jet=baseline.Jet(pt=30.0, eta=0.0, phi=0.0, mass=3.0, index=0),
            hard_fraction=1.0,
            hard_pt=30.0,
            total_pt=30.0,
            log_odds=-30.0,
        )
    ]
    recovered = baseline.append_graph_recovery_jets([], raw_records, args)
    assert recovered == []


# Check that latent graph recovery can add a high-raw-pT weak-posterior jet
def test_latent_graph_recovery_raw_graph_support() -> None:
    args = make_recovery_args(
        graph_kron_latent_recovery_raw_pt_min=70.0,
        graph_kron_latent_recovery_min_hard_pt=1.0,
        graph_kron_latent_recovery_max_add=1,
        graph_kron_latent_recovery_output_scale=0.30,
    )
    raw_records = [
        baseline.GraphJetRecord(
            jet=baseline.Jet(pt=80.0, eta=0.5, phi=-0.2, mass=8.0, index=0),
            hard_fraction=0.02,
            hard_pt=1.5,
            total_pt=80.0,
            log_odds=-5.0,
        )
    ]
    recovered = baseline.append_latent_graph_recovery_jets([], raw_records, args)
    assert len(recovered) == 1
    assert recovered[0].pt == pytest.approx(24.0)
    assert recovered[0].eta == pytest.approx(0.5)


# Check the shared massive-particle rapidity distance in the study adapter
def test_antikt_massive_particles_rapidity():
    particles = [baseline.WeightedParticle(20.0, 1.2, 0.0, 20.0),
                 baseline.WeightedParticle(5.0, 0.7, 0.0, 0.0)]
    jets = baseline.cluster_jets(particles, radius=0.4, pt_min=0.0)
    assert len(jets) == 1
    assert jets[0].pt == pytest.approx(25.0)
