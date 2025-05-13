import argparse
from collections import defaultdict

def parse_sam(sam_path, output_path):
    read_counts = defaultdict(int)

    # First pass: count read names
    with open(sam_path, "r") as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                continue  # skip header
            fields = line.strip().split("\t")
            read_id = fields[0]
            read_counts[read_id] += 1

    # Second pass: write read names and mapping type
    with open(sam_path, "r") as sam_file, open(output_path, "w") as out_file:
        seen = set()
        for line in sam_file:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            if read_id not in seen:
                tag = "unique" if read_counts[read_id] == 1 else "multi"
                out_file.write(f"{read_id}\t{tag}\n")
                seen.add(read_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract mapped read names and classify as unique or multi-mappers.")
    parser.add_argument("-i", "--input", required=True, help="Input SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output file to write read names and mapping status")

    args = parser.parse_args()
    parse_sam(args.input, args.output)

