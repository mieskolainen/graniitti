# Optimizer search state for distributed backends
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import random

import numpy as np

from core.tune import core as icetune_main
from core.tune.io import proposal_batch_size, proposals_independent
from core.tune.parameters.space import is_integer_bound, normalize_param_space, sample_config, typed_config


class SearchState:
    """Small in-memory optimizer rebuilt from completed and live trials."""

    # Initialize a search state for the configured typed parameter bounds
    def __init__(
        self,
        *,
        args,
        bounds: dict,
        initial_points: dict | None,
        parameter_topology: dict | None = None,
        async_proposals: bool = False,
    ):
        self.args = args
        self.bounds = normalize_param_space(bounds)
        self.initial_points = copy.deepcopy(initial_points)
        self.parameter_topology = copy.deepcopy(parameter_topology or {})
        self.async_proposals = bool(async_proposals)
        self.names = sorted(self.bounds)
        self.surrogate_fit_percentile = float(getattr(args, "surrogate_fit_percentile", 1.0))
        self._reset_optimizer_state()

    # Reset native optimizer state before replaying filtered trial history
    def _reset_optimizer_state(self) -> None:
        self.completed = set()
        self.completed_configs = {}
        self.completed_keys = set()
        self.observed = set()
        self.observed_keys = set()
        self.failed_observed = set()
        self.failed_configs = {}
        self.failed_keys = set()
        self.icebo = None
        self.amplitude = None
        self.amplitude_records = []
        self.amplitude_failures = []
        self.hebo = None
        self.hyperopt_base = None
        self.hyperopt_domain = None
        self.hyperopt_rng = None
        self.hyperopt_status_fail = None
        self.hyperopt_status_ok = None
        self.hyperopt_tpe = None
        self.hyperopt_trials = None
        self.hyperopt_live = {}
        self.optuna_study = None
        self.optuna_distributions = None
        self.optuna_live = {}
        self.optuna_pending = {}
        self.pending_configs = []
        self.pending_keys = set()
        self.icebo_batch_records = {}
        self.icebo_updated_batches = set()
        self.issued_proposals = set()
        self.issued_icebo_count = 0
        if self.args.algorithm == "ampfit" and self.bounds:
            from core.tune.optimizers.ampfit.optimizer import AmplitudeSearch

            self.amplitude = AmplitudeSearch(
                bounds=self.bounds, initial_points=self.initial_points,
                settings=self.args.ampfit_settings, seed=int(getattr(self.args, "rngseed", 0)))
        elif self.args.algorithm == "icebo" and self.bounds:
            self._init_icebo()
        elif self.args.algorithm == "hebo" and self.bounds:
            self._init_hebo()
        elif self.args.algorithm == "hyperopt" and self.bounds:
            self._init_hyperopt()
        elif self.args.algorithm == "optuna" and self.bounds:
            self._init_optuna()

    # Compute the completed records accepted for surrogate fitting
    def _surrogate_fit_records(self, records: list[dict]) -> list[dict]:
        percentile = 1.0 if self.args.algorithm == "icebo" else self.surrogate_fit_percentile
        return icetune_main.select_surrogate_fit_records(records, self.args.cost, percentile)

    # Initialize the standalone ICEBO optimizer state
    def _init_icebo(self) -> None:
        from core.tune.optimizers.icebo.optimizer import ICEBO

        self.icebo = ICEBO(
            self.bounds,
            settings=copy.deepcopy(self.args.icebo_settings),
            parameter_topology=self.parameter_topology,
            seed=int(getattr(self.args, "rngseed", 0)),
            warmup=0,
        )

    # Initialize the native HEBO optimizer state
    def _init_hebo(self) -> None:
        from hebo.design_space.design_space import DesignSpace

        from core.tune.optimizers.hebo.optimizer import create_hebo

        design = DesignSpace().parse(
            [
                {
                    "name": key,
                    "type": "int" if is_integer_bound(spec) else "num",
                    "lb": spec["lower"],
                    "ub": spec["upper"],
                }
                for key, spec in sorted(self.bounds.items())
            ]
        )
        self.hebo = create_hebo(
            design,
            parameter_topology=self.parameter_topology,
            rand_sample=0,
            scramble_seed=int(getattr(self.args, "rngseed", 0)),
            settings=copy.deepcopy(getattr(self.args, "hebo_settings", None)),
        )

    # Initialize the native Hyperopt optimizer state
    def _init_hyperopt(self) -> None:
        from hyperopt import STATUS_FAIL, STATUS_OK, Trials, base, hp, tpe

        self.hyperopt_base = base
        self.hyperopt_status_fail = STATUS_FAIL
        self.hyperopt_status_ok = STATUS_OK
        self.hyperopt_tpe = tpe
        self.hyperopt_rng = np.random.RandomState(int(getattr(self.args, "rngseed", 0)))
        space = {}
        for key, spec in sorted(self.bounds.items()):
            if is_integer_bound(spec):
                space[key] = hp.quniform(key, int(spec["lower"]), int(spec["upper"]), 1)
            else:
                space[key] = hp.uniform(key, spec["lower"], spec["upper"])
        self.hyperopt_domain = base.Domain(lambda spc: spc, space)
        self.hyperopt_trials = Trials()

    # Initialize the native Optuna optimizer state
    def _init_optuna(self) -> None:
        import optuna
        from optuna.distributions import FloatDistribution, IntDistribution

        sampler = optuna.samplers.TPESampler(
            seed=int(getattr(self.args, "rngseed", 0)),
            n_startup_trials=(0 if self.async_proposals else max(0, int(getattr(self.args, "rand_trials", 0)))),
            multivariate=False,
            group=False,
            constant_liar=True,
        )
        self.optuna_distributions = {}
        for key, spec in sorted(self.bounds.items()):
            if is_integer_bound(spec):
                self.optuna_distributions[key] = IntDistribution(int(spec["lower"]), int(spec["upper"]))
            else:
                self.optuna_distributions[key] = FloatDistribution(spec["lower"], spec["upper"])
        self.optuna_study = optuna.create_study(direction="minimize", sampler=sampler)

    # Compute a config coerced to the normalized parameter types
    def _typed_config(self, config: dict) -> dict:
        return typed_config(config, self.bounds)

    # Compute one hashable optimizer configuration identity
    def _config_key(self, config: dict) -> tuple:
        typed = self._typed_config(config)
        return tuple((name, typed[name]) for name in self.names)

    # Reconstruct optimizer-visible pending trials from scheduler state
    def set_pending(self, configs: list[dict]) -> None:
        self.pending_configs = [self._typed_config(config) for config in configs]
        self.pending_keys = {self._config_key(config) for config in self.pending_configs}
        if self.icebo is not None:
            self.icebo.set_pending(self.pending_configs)
            self.icebo.set_excluded(list(self.failed_configs.values()))
            return
        if self.optuna_study is None:
            return
        import optuna

        for typed in self.pending_configs:
            key = self._config_key(typed)
            live_keys = {self._config_key(config) for _, config in self.optuna_live.values()}
            if key in self.optuna_pending or key in live_keys:
                continue
            self.optuna_study.add_trial(
                optuna.trial.create_trial(
                    state=optuna.trial.TrialState.RUNNING, params=typed, distributions=self.optuna_distributions
                )
            )
            self.optuna_pending[key] = self.optuna_study.trials[-1].number

    # Restore issued proposal sequence state from the immutable scheduler journal
    def restore_issued(self, records: list[dict]) -> None:
        hyperopt_docs = []
        restore_icebo = self.icebo is not None and not self.issued_proposals
        latest_icebo = None
        for record in records:
            trial_id = str(record.get("trial_id") or "")
            payload = record.get("search_payload") or {}
            if not trial_id or not isinstance(payload, dict):
                continue
            proposal_index = payload.get("proposal_index")
            if proposal_index is not None:
                self.issued_icebo_count = max(self.issued_icebo_count, int(proposal_index) + 1)
            icebo_state = payload.get("icebo_state")
            if isinstance(icebo_state, dict):
                self.issued_icebo_count = max(self.issued_icebo_count, int(icebo_state.get("proposal_count", 0)))
                if latest_icebo is None or icebo_state["proposal_count"] > latest_icebo["proposal_count"]:
                    latest_icebo = icebo_state
            if trial_id in self.issued_proposals:
                continue
            if self.hyperopt_rng is not None and payload.get("kind") == "hyperopt":
                expected_seed = int(self.hyperopt_rng.randint(2**31 - 1))
                if expected_seed != int(payload.get("seed", -1)):
                    raise ValueError(f'Hyperopt proposal seed mismatch for trial "{trial_id}"')
                if trial_id not in self.completed:
                    doc = self._hyperopt_record_doc(record)
                    if doc is not None:
                        if trial_id in self.failed_observed:
                            doc["state"] = self.hyperopt_base.JOB_STATE_DONE
                            doc["result"] = {"status": self.hyperopt_status_fail}
                        hyperopt_docs.append(doc)
            if self.optuna_study is not None and trial_id not in self.completed:
                self._restore_optuna_record(record, failed=trial_id in self.failed_observed)
            self.issued_proposals.add(trial_id)
        if hyperopt_docs:
            self.hyperopt_trials.insert_trial_docs(hyperopt_docs)
            self.hyperopt_trials.refresh()
            self.hyperopt_live.update(
                {int(doc["tid"]): doc for doc in hyperopt_docs if doc["state"] != self.hyperopt_base.JOB_STATE_DONE}
            )
        if self.icebo is not None:
            if restore_icebo and latest_icebo is not None:
                self.icebo.load_state_dict(latest_icebo)
            self.icebo.set_replay_proposal_count(self.issued_icebo_count)

    # Generate a deterministic uniform random configuration
    def _uniform(self, index: int) -> dict:
        rng = random.Random(int(getattr(self.args, "rngseed", 0)) + int(index))
        return sample_config(bounds=self.bounds, uniform=rng.uniform, integer=rng.randint)

    # Compute the completed trial cost as a finite floating point value
    def _record_cost(self, record: dict) -> float:
        return float(record["metrics"][self.args.cost])

    # Compute a Hyperopt trial id for a completed tuning record
    def _hyperopt_tid(self, record: dict) -> int:
        payload = record["search_payload"]
        if payload.get("kind") == "hyperopt":
            return int(payload["tid"])
        return int(str(record["trial_id"]).rsplit("-", 1)[-1])

    # Build the Hyperopt misc payload for a known configuration
    def _hyperopt_misc(self, *, tid: int, config: dict) -> dict:
        return {
            "tid": tid,
            "cmd": ("domain_attachment", "FMinIter_Domain"),
            "workdir": None,
            "idxs": {key: [tid] for key in self.names},
            "vals": {key: [float(config[key])] for key in self.names},
        }

    # Rebuild one Hyperopt document from scheduler state
    def _hyperopt_record_doc(self, record: dict) -> dict | None:
        tid = self._hyperopt_tid(record)
        if any(int(trial["tid"]) == tid for trial in self.hyperopt_trials.trials):
            return None
        return self.hyperopt_trials.new_trial_docs(
            [tid], [None], [self.hyperopt_domain.new_result()], [self._hyperopt_misc(tid=tid, config=record["config"])]
        )[0]

    # Insert or update one completed Hyperopt trial from tuning history
    def _observe_hyperopt_record(self, record: dict) -> None:
        tid = self._hyperopt_tid(record)
        doc = self.hyperopt_live.pop(tid, None)
        if doc is None:
            doc = next((trial for trial in self.hyperopt_trials.trials if int(trial["tid"]) == tid), None)
        result = {"loss": self._record_cost(record), "status": self.hyperopt_status_ok}
        if doc is None:
            doc = self.hyperopt_trials.new_trial_docs(
                [tid], [None], [result], [self._hyperopt_misc(tid=tid, config=record["config"])]
            )[0]
            self.hyperopt_trials.insert_trial_docs([doc])
        doc["state"] = self.hyperopt_base.JOB_STATE_DONE
        doc["result"] = result
        self.hyperopt_trials.refresh()

    # Insert or update one completed Optuna trial from tuning history
    def _observe_optuna_record(self, record: dict) -> None:
        payload = record.get("search_payload", {})
        self._restore_optuna_rejected(payload)
        trial_index = payload.get("trial_index") if isinstance(payload, dict) else None
        value = self._record_cost(record)
        live = self.optuna_live.pop(int(trial_index), None) if trial_index is not None else None
        if live is not None:
            trial, proposed_config = live
            if self._typed_config(record["config"]) != proposed_config:
                raise ValueError(f'Optuna proposal identity mismatch for trial "{record.get("trial_id")}"')
            self.optuna_study.tell(trial, value)
            return

        pending_number = self.optuna_pending.pop(self._config_key(record["config"]), None)
        if pending_number is not None:
            self.optuna_study.tell(pending_number, value)
            return

        import optuna

        self.optuna_study.add_trial(
            optuna.trial.create_trial(
                params=self._typed_config(record["config"]), distributions=self.optuna_distributions, value=value
            )
        )

    # Rebuild one pending or failed Optuna trial from scheduler state
    def _restore_optuna_record(self, record: dict, *, failed: bool) -> None:
        import optuna

        self._restore_optuna_rejected(record.get("search_payload") or {})
        config = self._typed_config(record["config"])
        key = self._config_key(config)
        if key in self.optuna_pending:
            return
        state = optuna.trial.TrialState.FAIL if failed else optuna.trial.TrialState.RUNNING
        self.optuna_study.add_trial(
            optuna.trial.create_trial(state=state, params=config, distributions=self.optuna_distributions)
        )
        if not failed:
            self.optuna_pending[key] = self.optuna_study.trials[-1].number

    # Replay Optuna candidates rejected because they matched an existing config
    def _restore_optuna_rejected(self, payload: dict) -> None:
        if self.optuna_study is None:
            return
        import optuna

        for config in payload.get("rejected", []):
            self.optuna_study.add_trial(
                optuna.trial.create_trial(
                    state=optuna.trial.TrialState.FAIL,
                    params=self._typed_config(config),
                    distributions=self.optuna_distributions,
                )
            )

    # Observe completed trial records in trial-id order
    def observe(self, records: list[dict]) -> None:
        if self.amplitude is not None:
            self.amplitude_records = copy.deepcopy(records)
            return
        records = sorted(records, key=lambda record: str(record.get("trial_id") or ""))
        completed = {
            str(record["trial_id"])
            for record in records
            if isinstance(record, dict) and record.get("trial_id") is not None
        }
        fit_records = self._surrogate_fit_records(records)
        rebuild = self.args.algorithm != "icebo" and self.surrogate_fit_percentile < 1.0 and completed != self.completed
        if rebuild:
            self._reset_optimizer_state()
            fresh = fit_records
        else:
            fresh = [r for r in fit_records if r["trial_id"] not in self.observed]
        self.completed = completed
        self.completed_configs = {
            self._config_key(record["config"]): self._typed_config(record["config"]) for record in records
        }
        self.completed_keys = set(self.completed_configs)
        if not fresh:
            return
        if self.icebo is not None:
            self._observe_icebo_records(fresh)
        elif self.hebo is not None:
            import pandas as pd

            self.hebo.observe(
                pd.DataFrame([self._typed_config(r["config"]) for r in fresh]),
                np.array([[float(r["metrics"][self.args.cost])] for r in fresh]),
            )
        elif self.hyperopt_trials is not None:
            for record in fresh:
                self._observe_hyperopt_record(record)
        elif self.optuna_study is not None:
            for record in fresh:
                self._observe_optuna_record(record)
        for record in fresh:
            self.observed.add(record["trial_id"])
            key = self._config_key(record["config"])
            self.observed_keys.add(key)
            self.failed_keys.discard(key)
            self.failed_configs.pop(key, None)

    # Resolve failed native optimizer trials without fitting them
    def observe_failures(self, records: list[dict]) -> None:
        if self.amplitude is not None:
            self.amplitude_failures = copy.deepcopy(records)
            return
        fresh = [
            record
            for record in sorted(records, key=lambda item: str(item.get("trial_id") or ""))
            if record.get("trial_id") not in self.failed_observed
        ]
        for record in fresh:
            key = self._config_key(record["config"])
            if key not in self.completed_keys:
                self.failed_keys.add(key)
                self.failed_configs[key] = self._typed_config(record["config"])
            payload = record.get("search_payload") or {}
            if self.hyperopt_trials is not None and payload.get("kind") == "hyperopt":
                tid = self._hyperopt_tid(record)
                doc = self.hyperopt_live.pop(tid, None)
                if doc is not None:
                    doc["state"] = self.hyperopt_base.JOB_STATE_DONE
                    doc["result"] = {"status": self.hyperopt_status_fail}
                    self.hyperopt_trials.refresh()
            elif self.optuna_study is not None:
                import optuna

                trial_index = payload.get("trial_index") if payload.get("kind") == "optuna" else None
                live = self.optuna_live.pop(int(trial_index), None) if trial_index is not None else None
                if live is not None:
                    self.optuna_study.tell(live[0], state=optuna.trial.TrialState.FAIL)
                else:
                    number = self.optuna_pending.pop(self._config_key(record["config"]), None)
                    if number is not None:
                        self.optuna_study.tell(number, state=optuna.trial.TrialState.FAIL)
            elif self.icebo is not None and payload.get("kind") == "icebo":
                self._accumulate_icebo_turbo_trial(record["trial_id"], payload, None)
            self.failed_observed.add(record["trial_id"])

    # Replay ICEBO observations in their original proposal batches
    def _observe_icebo_records(self, records: list[dict]) -> None:
        grouped = {}
        for record in sorted(records, key=lambda item: str(item["trial_id"])):
            payload = record.get("search_payload", {})
            batch_id = payload.get("proposal_batch_id") if isinstance(payload, dict) else None
            key = ("batch", str(batch_id)) if batch_id is not None else ("trial", record["trial_id"])
            grouped.setdefault(key, []).append(record)
        proposal_count = 0
        for group in grouped.values():
            self.icebo.observe(
                [self._typed_config(record["config"]) for record in group],
                [float(record["metrics"][self.args.cost]) for record in group],
                [self._record_cost_error(record) for record in group],
                update_turbo=False,
            )
            for record in group:
                payload = record.get("search_payload", {})
                if not isinstance(payload, dict):
                    continue
                proposal_index = payload.get("proposal_index")
                if proposal_index is not None:
                    proposal_count = max(proposal_count, int(proposal_index) + 1)
                self._accumulate_icebo_turbo_trial(record["trial_id"], payload, self._typed_config(record["config"]))
        self.icebo.set_replay_proposal_count(proposal_count)

    # Update TuRBO only when every member of one proposal batch is complete
    def _accumulate_icebo_turbo_trial(self, trial_id: str, payload: dict, config: dict | None) -> None:
        if payload.get("kind") != "icebo" or payload.get("acquisition") == "sobol_warmup":
            return
        batch_id = payload.get("proposal_batch_id")
        batch_size = payload.get("proposal_batch_size")
        if batch_id is None or batch_size is None:
            if config is not None:
                self.icebo.update_turbo_batch([config], batch_id=f"trial:{trial_id}")
            return
        key = str(batch_id)
        if key in self.icebo_updated_batches:
            return
        expected = int(batch_size)
        if expected < 1:
            raise ValueError("ICEBO proposal batch size must be positive")
        batch = self.icebo_batch_records.setdefault(key, {})
        batch[str(trial_id)] = config
        if len(batch) < expected:
            return
        configs = [value for value in batch.values() if value is not None]
        if configs:
            self.icebo.update_turbo_batch(configs, batch_id=key)
        self.icebo_updated_batches.add(key)
        self.icebo_batch_records.pop(key, None)

    # Compute an optional scalar uncertainty for one optimizer observation
    def _record_cost_error(self, record: dict) -> float:
        metrics = record.get("metrics", {})
        for key in (f"{self.args.cost}_error", f"{self.args.cost}_uncertainty", "objective_error"):
            value = metrics.get(key)
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric) and numeric >= 0.0:
                return numeric
        likelihood = record.get("likelihood") or {}
        objective = likelihood.get("objective") or {}
        uncertainty = likelihood.get("mc_uncertainty") or {}
        if objective.get("name") == self.args.cost:
            value = uncertainty.get("sigma_objective")
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                numeric = -1.0
            if np.isfinite(numeric) and numeric >= 0.0:
                return numeric
        return 0.0

    # Ask Hyperopt for the next trial configuration
    def _ask_hyperopt(self, index: int) -> tuple[dict, dict]:
        seed = int(self.hyperopt_rng.randint(2**31 - 1))
        docs = self.hyperopt_tpe.suggest(
            [int(index)],
            self.hyperopt_domain,
            self.hyperopt_trials,
            seed,
            n_startup_jobs=(0 if self.async_proposals else max(0, int(getattr(self.args, "rand_trials", 0)))),
            verbose=False,
        )
        self.hyperopt_trials.insert_trial_docs(docs)
        self.hyperopt_trials.refresh()
        doc = docs[0]
        tid = int(doc["tid"])
        self.hyperopt_live[tid] = doc
        self.issued_proposals.add(f"trial-{int(index):06d}")
        return (
            self._typed_config({key: doc["misc"]["vals"][key][0] for key in self.names}),
            {"kind": "hyperopt", "seed": seed, "tid": tid},
        )

    # Ask Optuna for the next trial configuration
    def _ask_optuna(self, index: int) -> tuple[dict, dict]:
        trial = self.optuna_study.ask(fixed_distributions=self.optuna_distributions)
        config = self._typed_config(trial.params)
        self.optuna_live[int(index)] = (trial, config)
        self.issued_proposals.add(f"trial-{int(index):06d}")
        return (config, {"kind": "optuna", "number": int(trial.number), "trial_index": int(index)})

    # Ask Optuna again when its adaptive sampler repeats a blocked point
    def _ask_optuna_batch(self, index: int, count: int) -> list[tuple[dict, dict]]:
        blocked = self.completed_keys | self.failed_keys | self.pending_keys
        selected = set()
        rejected = []
        proposals = []
        limit = max(32, 8 * int(count))
        for _ in range(limit):
            trial_index = int(index) + len(proposals)
            config, payload = self._ask_optuna(trial_index)
            key = self._config_key(config)
            if key in blocked or key in selected:
                trial, _ = self.optuna_live.pop(trial_index)
                import optuna

                self.optuna_study.tell(trial, state=optuna.trial.TrialState.FAIL)
                rejected.append(config)
                continue
            if rejected:
                payload["rejected"] = copy.deepcopy(rejected)
                rejected.clear()
            proposals.append((config, payload))
            selected.add(key)
            if len(proposals) == int(count):
                return proposals
        raise RuntimeError("Optuna exhausted distinct adaptive configurations")

    # Compute one deterministic random proposal and its replay seed
    def _ask_random(self, index: int) -> tuple[dict, dict]:
        seed = int(getattr(self.args, "rngseed", 0)) + int(index)
        return self._uniform(index), {"kind": "random", "seed": seed}

    # Compute deterministic random proposals without invoking a surrogate fit
    def ask_random_many(self, index: int, count: int) -> list[tuple[dict, dict]]:
        target = max(0, int(count))
        proposals = []
        selected = set()
        attempt = 0
        while len(proposals) < target and attempt < 10000:
            config, payload = self._ask_random(int(index) + attempt)
            attempt += 1
            key = self._config_key(config)
            if key in self.completed_keys or key in self.pending_keys or key in self.failed_keys or key in selected:
                continue
            proposals.append((config, payload))
            selected.add(key)
        if len(proposals) != target:
            raise RuntimeError("random search exhausted distinct configurations")
        return proposals

    # Compute the configured constant lie without changing completed observations
    def _hebo_lie(self) -> float:
        settings = self.hebo.settings["pending"]
        values = self.hebo.y
        if settings["lie"] == "best":
            return float(np.min(values))
        if settings["lie"] == "mean":
            return float(np.mean(values))
        worst = float(np.max(values))
        penalty = max(float(np.std(values)), abs(worst) * 1.0e-6, 1.0e-6)
        return worst + settings["penalty_scale"] * penalty

    # Ask HEBO for one pending-aware adaptive block
    def _ask_hebo_batch(self, index: int, count: int) -> list[tuple[dict, dict]]:
        import pandas as pd

        proposals = []
        selected = set()
        excluded = {
            self._config_key(config): config for config in [*self.pending_configs, *self.failed_configs.values()]
        }
        excluded.update(
            {key: config for key, config in self.completed_configs.items() if key not in self.observed_keys}
        )
        optimizer = self.hebo
        synthetic_trials = 0
        if self.async_proposals and excluded and len(self.hebo.y):
            optimizer = copy.deepcopy(self.hebo)
            synthetic_trials = len(excluded)
            optimizer._model_config["synthetic_trials"] = synthetic_trials
            optimizer.observe(pd.DataFrame(list(excluded.values())), np.full((len(excluded), 1), self._hebo_lie()))
        base_seed = int(getattr(self.args, "rngseed", 0)) + int(index)
        for attempt in range(4):
            seed = base_seed + attempt
            optimizer.acquisition_seed = seed
            block = max(int(count) - len(proposals), 3)
            recommendations = optimizer.suggest(n_suggestions=block)
            accepted = []
            for offset, (_, row) in enumerate(recommendations.iterrows()):
                config = self._typed_config({key: row[key] for key in sorted(self.bounds)})
                key = self._config_key(config)
                if key in self.completed_keys or key in self.pending_keys or key in self.failed_keys or key in selected:
                    continue
                proposals.append((config, {"kind": "hebo", "batch_index": attempt * block + offset}))
                accepted.append(config)
                selected.add(key)
                if len(proposals) == count:
                    return proposals
            if accepted and len(self.hebo.y):
                optimizer = copy.deepcopy(optimizer)
                synthetic_trials += len(accepted)
                optimizer._model_config["synthetic_trials"] = synthetic_trials
                optimizer.observe(pd.DataFrame(accepted), np.full((len(accepted), 1), self._hebo_lie()))
        return proposals

    # Ask ICEBO in bounded pending-aware blocks after one surrogate fit
    def _ask_icebo_batch(self, count: int) -> list[tuple[dict, dict]]:
        proposals = []
        remaining = max(0, int(count))
        block = proposal_batch_size(self.args)
        while remaining:
            requested = min(block, remaining)
            batch_id = self.icebo.state_dict()["proposal_count"]
            configs, diagnostics = self.icebo.suggest(requested, return_diagnostics=True)
            if not configs:
                break
            state = self.icebo.state_dict()
            proposals.extend(
                (
                    self._typed_config(config),
                    {
                        "kind": "icebo",
                        "proposal_batch_id": batch_id,
                        "proposal_batch_size": len(configs),
                        "icebo_state": state,
                        **diagnostic,
                    },
                )
                for config, diagnostic in zip(configs, diagnostics, strict=True)
            )
            remaining -= len(configs)
        return proposals

    # Keep every emitted task distinct from terminal, pending and same-batch configurations
    def _unique_proposals(self, proposals: list[tuple[dict, dict]], index: int, count: int) -> list[tuple[dict, dict]]:
        blocked = self.completed_keys | self.failed_keys | self.pending_keys
        selected = set()
        unique = []
        for config, payload in proposals:
            typed = self._typed_config(config)
            key = self._config_key(typed)
            if key in blocked or key in selected:
                continue
            unique.append((typed, payload))
            selected.add(key)
            if len(unique) == count:
                return unique
        random_only = bool(proposals) and all(
            (payload or {}).get("kind") in {"basic", "cold", "initial", "random"} for _, payload in proposals
        )
        if not random_only:
            raise RuntimeError(f"{self.args.algorithm} returned duplicate or blocked adaptive proposals")
        attempt = 0
        while len(unique) < count and attempt < 10000:
            fallback_index = 1_000_000 + int(index) + attempt
            config, payload = self._ask_random(fallback_index)
            attempt += 1
            key = self._config_key(config)
            if key in blocked or key in selected:
                continue
            unique.append((config, payload))
            selected.add(key)
        if len(unique) != count:
            raise RuntimeError("optimizer exhausted distinct configurations")
        return unique

    # Ask the active optimizer for an aligned batch of configurations and metadata
    def ask_many(self, index: int, count: int) -> list[tuple[dict, dict]]:
        if self.amplitude is not None:
            return self.amplitude.ask(records=self.amplitude_records, failures=self.amplitude_failures,
                                      pending=self.pending_configs, count=count, cost=self.args.cost)
        target = max(0, int(count))
        remaining = target
        if remaining == 0 or (self.pending_keys and not self.async_proposals):
            return []

        proposals = []
        next_index = int(index)
        if next_index == 0 and self.initial_points is not None:
            proposals.append(self.ask(next_index))
            next_index += 1
            remaining -= 1
        if remaining == 0:
            return self._unique_proposals(proposals, index, target)

        random_phase = proposals_independent(self.args, next_index)
        if self.icebo is not None:
            warmup = max(0, int(getattr(self.args, "rand_trials", 0)))
            random_count = min(remaining, max(0, warmup - next_index))
            proposals.extend(self.ask_random_many(next_index, random_count))
            next_index += random_count
            remaining -= random_count
            if remaining > 0:
                proposals.extend(self._ask_icebo_batch(remaining))
            if not proposals:
                return []
            target = len(proposals)
        elif self.hebo is not None:
            warmup = max(0, int(getattr(self.args, "rand_trials", 0)))
            random_count = min(remaining, max(0, warmup - next_index))
            proposals.extend(self.ask_random_many(next_index, random_count))
            next_index += random_count
            remaining -= random_count
            if remaining > 0:
                proposals.extend(self._ask_hebo_batch(next_index, remaining))
        elif self.async_proposals and random_phase:
            proposals.extend(self.ask_random_many(next_index, remaining))
        elif self.optuna_study is not None:
            proposals.extend(self._ask_optuna_batch(next_index, remaining))
        else:
            proposals.extend(self.ask(next_index + offset) for offset in range(remaining))
        return self._unique_proposals(proposals, index, target)

    # Ask the optimizer for the next trial configuration
    def ask(self, index: int) -> tuple[dict, dict]:
        if index == 0 and self.initial_points is not None:
            config = self._typed_config(copy.deepcopy(self.initial_points))
            if self.icebo is not None:
                self.icebo.reserve(config)
            return config, {"kind": "initial"}
        if self.hyperopt_trials is not None:
            return self._ask_hyperopt(index)
        if self.optuna_study is not None:
            return self._ask_optuna(index)
        if self.icebo is not None and index < int(getattr(self.args, "rand_trials", 0)):
            return self._ask_random(index)
        if self.icebo is not None:
            proposals = self._ask_icebo_batch(1)
            if not proposals:
                raise RuntimeError("ICEBO is waiting for pending evaluations, use ask_many to poll")
            return proposals[0]
        if self.hebo is not None and index >= int(getattr(self.args, "rand_trials", 0)):
            proposals = self._ask_hebo_batch(index, 1)
            if proposals:
                return proposals[0]
            raise RuntimeError("HEBO returned no distinct proposal")
        return self.ask_random_many(index, 1)[0]
