import os
import glob
import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run MARIA on chunked files (serial execution)."
    )

    parser.add_argument(
        "--split_folder",
        required=True,
        help="Folder containing the chunked input files."
    )

    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix used in chunk filenames (e.g., NoLinker_V1)."
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        required=True,
        help="Percentile cutoff for MARIA (e.g., 63)."
    )

    parser.add_argument(
        "--tpm_reference",
        default=None,
        help="Optional TPM reference (e.g., K562)."
    )

    return parser.parse_args()


def run_maria_on_chunk(chunk_file, cutoff, tpm_reference):
    """Run MARIA on a single chunk file (serial execution)."""
    print("Running MARIA on {}".format(chunk_file))

    cmd = ["python", "maria.py", chunk_file, "-cut_off", str(cutoff)]

    if tpm_reference is not None:
        cmd += ["-tpm_reference", tpm_reference]

    subprocess.call(cmd)


def main():
    args = parse_args()

    split_folder = args.split_folder
    output_prefix = args.output_prefix
    cutoff = args.cutoff
    tpm_reference = args.tpm_reference

    # Find chunk files like:  <split_folder>/<output_prefix>_chunk_1.txt
    pattern = os.path.join(split_folder, "{}_chunk_*.txt".format(output_prefix))
    chunk_files = sorted(glob.glob(pattern))

    if not chunk_files:
        print("No chunk files found with pattern: {}".format(pattern))
        return

    print("Found {} chunk files.".format(len(chunk_files)))

    # Serial execution
    for chunk_file in chunk_files:
        run_maria_on_chunk(chunk_file, cutoff, tpm_reference)

    print("All chunks processed.")


if __name__ == "__main__":
    main()