import os
import random
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path


def load_dataset(dataset_path: Path) -> list[tuple[int, str, str]]:
    lines = dataset_path.read_text().splitlines()
    dataset = []
    for i in range(len(lines) // 3):
        dataset.append((int(lines[3 * i]), lines[3 * i + 1], lines[3 * i + 2]))
    return dataset


def make_estimator_process(estimator_path: Path) -> subprocess.Popen:
    return subprocess.Popen(
        ["uv", "run", estimator_path.as_posix()],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def estimate_mutations(
    estimator_path: Path, original_dna: str, mutated_dna: str
) -> int:
    with make_estimator_process(estimator_path) as proc:
        proc.stdin.write(f"{original_dna}\n")
        proc.stdin.write(f"{mutated_dna}\n")
        estimate = int(proc.stdout.readline())
    return estimate


def mse(xs: list[int], ys: list[int]) -> float:
    n = len(xs)
    return (1 / n) * sum((x - y) ** 2 for x, y in zip(xs, ys))


def mean(xs: list[int]) -> float:
    return sum(xs) / len(xs)


def std(xs: list[int]) -> float:
    mu = mean(xs)
    return (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5


def run_estimator(
    estimator_path: Path,
    dataset: list[tuple[int, str, str]],
) -> list[int]:
    cpus = os.cpu_count()
    with Pool(cpus - 1) as pool:
        estimates = pool.starmap(
            estimate_mutations,
            [
                (estimator_path, original_dna, mutated_dna)
                for _, original_dna, mutated_dna in dataset
            ],
        )

    return estimates


def print_statistics(dataset: list[tuple[int, str, str]], estimates: list[int]) -> None:
    random.seed(0)

    actuals = [m for m, _, _ in dataset]

    print("OVERALL".center(80, "*"))
    print()
    print(f"MSE: {mse(actuals, estimates):.2f}")
    print()

    actual_to_estimate = {}
    for actual, estimate in zip(actuals, estimates):
        if actual not in actual_to_estimate:
            actual_to_estimate[actual] = []
        actual_to_estimate[actual].append(estimate)

    print("INDIVIDUAL".center(80, "*"))
    print()
    print("Mutations     : MSE     Mean    Std     Min     Max")
    for actual, _estimates in sorted(actual_to_estimate.items()):
        n = len(_estimates)
        _actuals = [actual for _ in range(n)]
        print(
            f"{actual:<14}: {mse(_actuals, _estimates):<7.2f} {mean(_estimates):<7.2f} {std(_estimates):<7.2f} {min(_estimates):<7} {max(_estimates):<7}"
        )

    print()
    print("Mutations\t - Actual number of mutations")
    print("MSE\t\t - Mean squared error")
    print("Mean\t\t - Mean estimate")
    print("Std\t\t - Standard deviation of estimates")
    print("Min\t\t - Minimum estimate")
    print("Max\t\t - Maximum estimate")
    print()

    print("GOOD ESTIMATES".center(80, "*"))
    print()
    min_error = float("inf")
    best_estimates = []
    for (actual, original, mutated), estimate in zip(dataset, estimates):
        error = abs(actual - estimate)
        if error > min_error:
            continue
        if error < min_error:
            min_error = error
            best_estimates.clear()
        best_estimates.append((actual, estimate, original, mutated))

    random.shuffle(best_estimates)
    best_estimates = best_estimates[:10]

    for actual, estimate, original, mutated in best_estimates:
        print(f"Mutations: {actual}")
        print(f"Estimate : {estimate}")
        print(original)
        print(mutated)
        print()

    print("BAD ESTIMATES".center(80, "*"))
    print()

    max_error = float("-inf")
    worst_estimates = []
    for (actual, original, mutated), estimate in zip(dataset, estimates):
        error = abs(actual - estimate)
        if error < max_error:
            continue
        if error > max_error:
            max_error = error
            worst_estimates.clear()
        worst_estimates.append((actual, estimate, original, mutated))

    random.shuffle(worst_estimates)
    worst_estimates = worst_estimates[:10]

    for actual, estimate, original, mutated in worst_estimates:
        print(f"Mutations: {actual}")
        print(f"Estimate : {estimate}")
        print(original)
        print(mutated)
        print()

    print("RANDOM ESTIMATES".center(80, "*"))
    print()

    random_estimates = [
        (actual, estimate, original, mutated)
        for (actual, original, mutated), estimate in zip(dataset, estimates)
    ]
    random.shuffle(random_estimates)
    random_estimates = random_estimates[:10]

    for actual, estimate, original, mutated in random_estimates:
        print(f"Mutations: {actual}")
        print(f"Estimate : {estimate}")
        print(original)
        print(mutated)
        print()


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: uv run evaluate.py <dataset_path> <estimator_path>")
        return

    dataset_path = Path(sys.argv[1])
    estimator_path = Path(sys.argv[2])

    if dataset_path.suffix == ".py" or estimator_path.suffix == ".txt":
        print("Usage: uv run evaluate.py <dataset_path> <estimator_path>")
        return

    dataset = load_dataset(dataset_path)

    estimates = run_estimator(estimator_path, dataset)
    print_statistics(dataset, estimates)


if __name__ == "__main__":
    main()
