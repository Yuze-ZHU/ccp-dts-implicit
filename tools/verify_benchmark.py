import argparse
from pathlib import Path


def numeric_rows(path, minimum_columns):
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < minimum_columns:
            continue
        try:
            rows.append([float(value) for value in fields])
        except ValueError:
            continue
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--result",
        type=Path,
        default=Path("results/lymberopoulos1993/NeNiEeAveDT.dat"),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("reference/benchmark_values.dat"),
    )
    parser.add_argument("--tolerance", type=float, default=5.0e-3)
    args = parser.parse_args()

    result_rows = numeric_rows(args.result, 8)
    reference_rows = numeric_rows(args.reference, 3)
    if not result_rows:
        raise SystemExit(f"No numerical rows found in {args.result}")
    if not reference_rows:
        raise SystemExit(f"No numerical rows found in {args.reference}")

    observed = result_rows[-1][5:8]
    expected = reference_rows[-1][:3]
    names = ("NeAve", "NiAve", "EeAve")
    errors = [abs(value - target) / abs(target) for value, target in zip(observed, expected)]

    for name, value, target, error in zip(names, observed, expected, errors):
        print(
            f"{name}: observed={value:.8e} expected={target:.8e} "
            f"relative_error={error:.6e}"
        )

    maximum = max(errors)
    print(f"maximum_relative_error={maximum:.6e}")
    if maximum > args.tolerance:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
