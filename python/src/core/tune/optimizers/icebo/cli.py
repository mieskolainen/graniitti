#!/usr/bin/env python3
# Run reproducible ICEBO optimizer and scaling benchmarks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import json

from core.tune.optimizers.icebo.benchmark import benchmark_optimizers, benchmark_scaling, save_benchmark


# Parse one comma-separated integer sequence
def integer_list(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


# Parse one comma-separated string sequence
def string_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


# Run paired optimizer and optional scaling benchmarks
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark ICEBO on paired deterministic and noisy problems"
    )
    parser.add_argument("--optimizers", type=string_list, default=["icebo", "hebo", "optuna"])
    parser.add_argument(
        "--problems",
        type=string_list,
        help="Comma-separated benchmark problem names; defaults to the complete suite",
    )
    parser.add_argument("--seeds", type=integer_list, default=[11, 23, 47])
    parser.add_argument("--budget", type=int, default=30)
    parser.add_argument("--initial_points", type=int)
    parser.add_argument("--batch_size", type=int, default=1)
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--scaling_sizes", type=integer_list, default=[])
    parser.add_argument("--scaling_dimension", type=int, default=12)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--output", default="tmp/icebo_benchmark.json")
    args = parser.parse_args()

    result = benchmark_optimizers(
        optimizers=args.optimizers,
        seeds=args.seeds,
        budget=args.budget,
        initial_points=args.initial_points,
        fast=args.fast,
        problem_names=args.problems,
        batch_size=args.batch_size,
    )
    if args.scaling_sizes:
        result["scaling"] = benchmark_scaling(
            sizes=args.scaling_sizes,
            dimension=args.scaling_dimension,
            seed=args.seeds[0],
            device=args.device,
        )
    save_benchmark(result, args.output)
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    print(f"Saved benchmark: {args.output}")


if __name__ == "__main__":
    main()
