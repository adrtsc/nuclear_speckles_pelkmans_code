from pathlib import Path
import pandas as pd
import re
import numpy as np
from Bio.SeqUtils import GC
from tqdm import tqdm
from Bio import SeqIO
import sys


def tx_info(x):
    g_name = x.attribute.split("gene_id ")[1].split(";")[0].replace('"', '')
    g_symbol = x.attribute.split("gene_name ")[1].split(";")[0].replace('"', '')
    g_type = x.attribute.split("gene_type ")[1].split(";")[0].replace('"', '')
    if x.feature != "gene":
        tx_name = x.attribute.split(
            "transcript_id ")[1].split(";")[0].replace('"', '')
        tx_type = x.attribute.split(
            "transcript_type ")[1].split(";")[0].replace('"', '')
    else:
        tx_name = np.nan
        tx_type = np.nan
    return g_name, tx_name, g_type, tx_type, g_symbol


def tx_info_gff3(x):
    g_name = x.attribute.split("gene_id=")[1].split(";")[0].replace('"', '')
    g_symbol = x.attribute.split("gene_name=")[1].split(";")[0].replace('"', '')
    g_type = x.attribute.split("gene_type=")[1].split(";")[0].replace('"', '')
    if x.feature != "gene":
        tx_name = x.attribute.split(
            "transcript_id=")[1].split(";")[0].replace('"', '')
        tx_type = x.attribute.split(
            "transcript_type=")[1].split(";")[0].replace('"', '')
    else:
        tx_name = np.nan
        tx_type = np.nan
    return g_name, tx_name, g_type, tx_type, g_symbol


def get_tx_sequences(annotations, record_dict):
    transcript_list = []

    for tx_name in tqdm(annotations['tx_name'].dropna().unique(),
                        total=annotations['tx_name'].unique().shape[0]):
        sequences = []
        features = []
        # get all information about the transcript
        tx = annotations.loc[annotations.tx_name == tx_name]

        # define start of tx:
        tx_start = tx.loc[tx.feature == 'transcript'].start.values[0]

        # iterate over features and append name & sequence
        for ii, f in tx.iterrows():
            # append the feature name and feature sequence
            features.append(f['feature'])
            if f.strand == "+":
                f_seq = record_dict[f['seqname']][f.start:f.end].seq
                if f.feature != 'gene':
                    f_seq = f_seq.transcribe()
                else:
                    pass
            elif f.strand == "-":
                f_seq = record_dict[f['seqname']][f.start:f.end].seq
                if f.feature != 'gene':
                    f_seq = f_seq.reverse_complement().transcribe()
                else:
                    f_seq = f_seq.reverse_complement()

            sequences.append(str(f_seq))

        # second iteration to define intron sequences
        sequence_list = list(sequences[features == 'transcript'])

        for ii, f in tx.iterrows():
            # set all remaining features of transcript to N (introns will be left)
            if f.feature not in ['gene', 'transcript']:
                feature_start_id = f.start - tx_start
                feature_length = f.end - f.start
                feature_end_id = feature_start_id + feature_length

                sequence_list[
                feature_start_id:feature_end_id] = "N" * feature_length

        sequence = ''.join(map(str, sequence_list))

        # find all the introns from the masked sequence
        introns = re.findall(r'(?:(?<=N)|^)[^,N]+', sequence)
        sequences = sequences + introns
        features = features + ['intron' for intron in introns]


        # make a dataframe containing all the relevant sequences
        res = pd.DataFrame({'feature': features,
                            'sequence': sequences,
                            'tx_name': tx_name})

        # spliced sequence
        full_seq = res.loc[res.feature == 'transcript']['sequence'].values[0]

        for i, f in res.iterrows():
            if f.feature == 'intron':
                full_seq = full_seq.replace(f['sequence'], "")

        sequences = sequences + [full_seq]
        features = features + ['spliced_transcript']

        # make a dataframe containing all the relevant sequences
        res = pd.DataFrame({'feature': features,
                            'sequence': sequences,
                            'tx_name': tx_name})

        transcript_list.append(res)

    res = pd.concat(transcript_list)

    return res

def get_intron_exon_features(res):

    exon_intron = res.loc[res.feature.isin(['exon', 'intron'])]
    exon_intron_counts = exon_intron.groupby(['tx_name', 'feature']).aggregate(
        {'feature': 'count'}).unstack(fill_value=0)

    exon_intron_agg = exon_intron.groupby(['tx_name', 'feature']).aggregate(
        {'GC_content': ['min', 'max', 'std', 'mean'],
         'length': ['min', 'max', 'std', 'mean', 'sum'],
         'AGCCC_count': ['min', 'max', 'std', 'mean', 'sum'],
         'AGCCC_density': ['min', 'max', 'std', 'mean'],
         #'splice_site_count': ['min', 'max', 'std', 'mean', 'sum'],
         #'splice_site_density': ['min', 'max', 'std', 'mean'],
         'm6A_count': ['min', 'max', np.std, 'mean', 'sum'],
         'm6A_density': ['min', 'max', np.std, 'mean'],
         'GGUCCGUCAU_count': ['min', 'max', np.std, 'mean', 'sum'],
         'GGUCCGUCAU_density': ['min', 'max', np.std, 'mean']
         }).unstack()



    exon_intron_agg.columns = ['_'.join(col).strip() for col in
                               exon_intron_agg.columns.values]

    exon_intron_counts.columns = ['_count_'.join(col).strip() for col in
                             exon_intron_counts.columns.values]

    return exon_intron_agg,  exon_intron_counts


def get_UTR_features(res):

    UTR = res.loc[res.feature.isin(['five_prime_UTR', 'three_prime_UTR'])]

    UTR_agg = UTR.groupby(['tx_name', 'feature']).aggregate(
        {'GC_content': ['mean'],
         'length': ['sum'],
         'AGCCC_count': ['sum'],
         'AGCCC_density': ['mean'],
         # 'splice_site_count': ['min', 'max', 'std', 'mean', 'sum'],
         # 'splice_site_density': ['min', 'max', 'std', 'mean'],
         'm6A_count': ['sum'],
         'm6A_density': ['mean'],
         'GGUCCGUCAU_count': ['sum'],
         'GGUCCGUCAU_density': ['sum']
         }).unstack()

    UTR_agg.columns = ['_'.join(col).strip() for col in UTR_agg.columns.values]
    UTR_agg.columns = [col.replace("_sum", '') for col in
                       UTR_agg.columns.values]
    UTR_agg.columns = [col.replace("_mean", '') for col in
                       UTR_agg.columns.values]

    # density needs to be re-computed across all UTR regions (some tx have it split)
    UTR_agg['m6A_density_three_prime_UTR'] = UTR_agg['m6A_count_three_prime_UTR']/UTR_agg['length_three_prime_UTR']
    UTR_agg['m6A_density_five_prime_UTR'] = UTR_agg['m6A_count_five_prime_UTR']/UTR_agg['length_five_prime_UTR']
    UTR_agg['AGCCC_density_three_prime_UTR'] = UTR_agg['AGCCC_count_three_prime_UTR']/UTR_agg['length_three_prime_UTR']
    UTR_agg['AGCCC_density_five_prime_UTR'] = UTR_agg['AGCCC_count_five_prime_UTR']/UTR_agg['length_five_prime_UTR']
    UTR_agg['GGUCCGUCAU_density_three_prime_UTR'] = UTR_agg['GGUCCGUCAU_count_three_prime_UTR']/UTR_agg['length_three_prime_UTR']
    UTR_agg['GGUCCGUCAU_density_five_prime_UTR'] = UTR_agg[ 'GGUCCGUCAU_count_five_prime_UTR']/UTR_agg['length_five_prime_UTR']

    return UTR_agg


def get_tx_features(annotations, record_dict):

    res = get_tx_sequences(annotations, record_dict)
    
    regions_to_keep = ['transcript',
                       'exon',
                       'intron',
                       'five_prime_UTR',
                       'three_prime_UTR', 'spliced_transcript']
    
    res = res.loc[res.feature.isin(regions_to_keep)]

    # add features to result table
    res['GC_content'] = res['sequence'].transform(lambda x: GC(x))
    res['length'] = res['sequence'].transform(lambda x: len(x))

    # drop features that have length 0
    res = res.loc[res.length > 0]

    #res['AGCCC_count'] = res.apply(lambda x: x['sequence'].count("AGCCC"),
    #                               axis=1)
    #res['AGCCC_density'] = res.apply(lambda x: x['sequence'].count("AGCCC")/len(x['sequence']),
    #                               axis=1)
    res['AGCCC_count'] = res['sequence'].transform(lambda x: len(re.findall('[U, A].{4}[G, C].{2}AGCCC', str(x))))
    res['AGCCC_density'] = res['sequence'].transform(
        lambda x: len(re.findall('[U, A].{4}[G, C].{2}AGCCC', str(x))) / len(x))

    res['m6A_count'] = res['sequence'].transform(lambda x: x.count("GGACU"))
    res['m6A_density'] = res['sequence'].transform(lambda x: x.count("GGACU")/len(x))
    res['GGUCCGUCAU_count'] = res['sequence'].transform(lambda x: x.count("GGUCCGUCAU"))
    res['GGUCCGUCAU_density'] = res['sequence'].transform(lambda x: x.count("GGUCCGUCAU") / len(x))




    # drop duplicates
    res = res.drop_duplicates()

    # aggregate the results
    tx_features = res.loc[res.feature.isin(['transcript'])]

    exon_intron_agg, exon_intron_counts = get_intron_exon_features(res)

    UTR_agg = get_UTR_features(res)

    spliced_agg = res.loc[res.feature.isin(['spliced_transcript'])].set_index('tx_name')[['GC_content', 'length',
                                                                     'AGCCC_count', 'AGCCC_density',
                                                                     'm6A_count', 'm6A_density',
                                                                     'GGUCCGUCAU_count', 'GGUCCGUCAU_density']]
    spliced_agg.columns = [x + "_spliced_transcript" for x in spliced_agg.columns]

    all_features = pd.concat(
        [tx_features.set_index("tx_name"), exon_intron_counts,
         exon_intron_agg, UTR_agg, spliced_agg], axis=1)

    annot_agg = annotations.groupby(['gene_id', 'tx_name']).agg({
        'gene_type': lambda x: np.unique(x)[0],
        'tx_type': lambda x: np.unique(x)[0],
        'gene_name': lambda x: np.unique(x)[0]})

    annot_agg.reset_index(inplace=True)

    all_features = pd.merge(all_features, annot_agg, on='tx_name')

    all_features = all_features.drop(['feature', 'sequence'], axis=1)

    return all_features


def main():

    n_jobs = int(sys.argv[1])
    job_id = int(sys.argv[2])

    genome_path = Path(
        r"/data/active/atschan/20220714_apex2_sequencing/reference/"
        r"GRCh38.primary_assembly.genome.fa")

    #genome_path = Path(
    #    r"C:\Users\Adria\Desktop\work\reference\GRCh38.p14.genome.fa")

    annotation_path = Path(
        r"/data/active/atschan/20220714_apex2_sequencing/reference/"
        r"gencode.v43.primary_assembly.annotation.gff3")

     annotation_path = Path(r"Z:\20220714_apex2_sequencing\reference\gencode.v43.primary_assembly.annotation.gff3")

    # read the genome sequences
    print("reading genome")
    record_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    # read the annotations
    print("reading annotations")
    annotations = pd.read_table(annotation_path,
                                comment="#", sep="\t",
                                names=['seqname', 'source', 'feature',
                                       'start', 'end', 'score', 'strand',
                                       'frame', 'attribute'])

    # add attributes to df
    annotations['gene_id'], \
    annotations['tx_name'], \
    annotations['gene_type'], \
    annotations['tx_type'], \
    annotations['gene_name'] = zip(
        *annotations.apply(lambda x: tx_info_gff3(x), axis=1))

    all_tx = annotations.tx_name.unique()

    chunks = np.array_split(all_tx, n_jobs)
    data = annotations.loc[annotations.tx_name.isin(chunks[job_id])]

    t = get_tx_features(data, record_dict)

    t.to_csv(fr"/data/active/atschan/20220714_apex2_sequencing/fishpond_results/features/gff3/20230724/features_{job_id}.csv")


if __name__ == "__main__":
    main()
