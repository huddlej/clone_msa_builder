#!/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()

    new_records = []
    for record in SeqIO.parse(args.input_file, "fasta"):
        new_records.append(SeqRecord(Seq(str(record.seq)), id=record.id, description=""))

    # Write out all trimmed records to the given output file.
    SeqIO.write(new_records, args.output_file, "fasta")
