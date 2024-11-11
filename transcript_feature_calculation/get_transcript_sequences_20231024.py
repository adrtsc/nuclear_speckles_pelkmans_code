from pathlib import Path
import pandas as pd
from Bio import SeqIO


def main():


    sequence_path = Path(
        r"Z:\20220714_apex2_sequencing\reference\gencode.v44.pc_transcripts.fa")

    record_dict = SeqIO.to_dict(SeqIO.parse(sequence_path, "fasta"))

    names = []
    sequences = []

    for key, value in record_dict.items():
        sequences.append(str(value.seq))
        names.append(key.split("|")[0])

    data = pd.DataFrame({'tx_id': names, 'sequence': sequences})

    data['length'] = data.sequence.transform(lambda x: len(x))

    data.to_csv(Path(r"Z:\20220714_apex2_sequencing\transcript_sequences\protein_coding_transcripts_sequences.csv"))






if __name__ == "__main__":
    main()