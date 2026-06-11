import os
import argparse
import sys

def main():
    # ---- Parse command-line arguments ----
    parser = argparse.ArgumentParser(
        description="Combine one or more MARIA input files and split into chunks."
    )
    parser.add_argument(
        "input_filenames",
        nargs="+",
        help="One or more MARIA input text files to combine and split."
    )
    parser.add_argument(
        "output_prefix",
        help="Base name used for output files and split folder (split_<output_prefix>)."
    )
    parser.add_argument(
        "--max_samples_per_chunk",
        type=int,
        default=500000,
        help="Maximum number of data lines per chunk file (default: 500000)."
    )

    args = parser.parse_args()

    input_filenames = args.input_filenames
    output_prefix = args.output_prefix
    max_samples_per_chunk = args.max_samples_per_chunk

    # ---- Load and validate inputs ----
    combined_data = []
    header = None

    for fname in input_filenames:
        if not os.path.exists(fname):
            print("Input file not found: {}".format(fname))
            sys.exit(1)

        with open(fname, "r") as f:
            lines = f.readlines()

        if not lines:
            print("Input file is empty: {}".format(fname))
            sys.exit(1)

        file_header = lines[0]

        if header is None:
            header = file_header
        elif file_header != header:
            print("❌ Header mismatch detected in file:", fname)
            print("Expected header:")
            print(header.strip())
            print("Found header:")
            print(file_header.strip())
            sys.exit(1)

        combined_data.extend(lines[1:])  # append data lines only

    if not combined_data:
        print("No data rows found after combining input files.")
        sys.exit(1)

    # ---- Create output directory ----
    out_dir = "split_{}".format(output_prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # ---- Split combined data into chunks ----
    num_chunks = (len(combined_data) + max_samples_per_chunk - 1) // max_samples_per_chunk

    for i in range(num_chunks):
        chunk = combined_data[
            i * max_samples_per_chunk : (i + 1) * max_samples_per_chunk
        ]
        chunk_filename = os.path.join(
            out_dir, "{}_chunk_{}.txt".format(output_prefix, i + 1)
        )

        with open(chunk_filename, "w") as out_f:
            out_f.write(header)
            out_f.writelines(chunk)

    print(
        "Splitting completed. {} chunks created in '{}' from {} input files.".format(
            num_chunks, out_dir, len(input_filenames)
        )
    )

if __name__ == "__main__":
    main()