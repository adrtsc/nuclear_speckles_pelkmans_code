from pathlib import Path
import pandas as pd
import numpy as np

def load_dge_results(dge_path):
    # load dge results:
    dge_results = pd.read_csv(dge_path)
    cols_to_drop = [col for col in dge_results.columns if (".y" in col)]
    dge_results.drop(cols_to_drop, axis=1, inplace=True)
    dge_results.columns = [col.split('.x')[0] for col in dge_results.columns]

    dge_results['ensembl_transcript_id'] = dge_results['tx_name'].transform(
        lambda x: x.split(".")[0])

    return dge_results


def load_rate_data(rate_path):
    # load rate data
    rates_data = pd.read_excel(rate_path, sheet_name=[1, 2])
    rates = []
    for rep in rates_data:
        data = rates_data[rep]
        data['replicate'] = rep
        rates.append(data)

    rates = pd.concat(rates)
    rates = rates.groupby('Symbol').mean()

    # remove some rates (defined by previous exploratory analysis)
    rates_to_use = [c for c in rates.columns if
                    (('k_' in c) | ('Half_life' in c)) & ~('quantile' in c)]
    rates = rates[rates_to_use]
    rates = rates.drop(['k_nucexp_from_chr.Mean', 'Half_life_nucdeg.Mean',
                        'k_nucexp_from_chr.MAP',
                        'Half_life_nucexp_from_nucdeg.MAP',
                        'Half_life_poly_entry.MAP'], axis=1)

    return rates


def load_RBP_features(RBP_path):
    # load RBP features
    RBP_files = list(RBP_path.glob("*.csv"))
    RBP_features = [pd.read_csv(fyle) for fyle in RBP_files]
    RBP_features = pd.concat(RBP_features, axis=0)

    return RBP_features



def load_sequence_features(gff_path):
    # load sequence features
    feature_files = list(gff_path.glob("*.csv"))
    feature_list = [pd.read_csv(fyle) for fyle in feature_files]
    gff_features = pd.concat(feature_list)

    return gff_features


def load_TPM(TPM_path):
    # load salmon TPM
    TPM_data = pd.read_csv(TPM_path)
    TPM_data['tx_id'] = TPM_data['tx_id'].transform(lambda x: x.split(".")[0])

    return TPM_data[['tx_id', 'TPM']]


def load_promoter_motifs(promoter_path):
    # load promoter motifs
    promoter_motif_files = list(promoter_path.glob("*.fps"))
    promoter_motifs = []

    for fyle in promoter_motif_files:
        motif_name = fyle.stem.split("_")[1]
        motifs = pd.read_table(Path(fyle), header=None, delim_whitespace=True).rename(columns={1: 'SYMBOL'})
        motifs['SYMBOL'] = motifs['SYMBOL'].apply(lambda x: x.split("_")[0])

        counts = motifs.groupby('SYMBOL').size().rename(f'{motif_name}_count')
        promoter_motifs.append(counts)

    promoter_motifs = pd.concat(promoter_motifs, axis=1).fillna(0)

    # rename to density (will be changed to density below)
    promoter_motifs.columns = [c.replace("count", "density") for c in promoter_motifs.columns]

    return promoter_motifs


def main():

    # define output path of merged features
    output_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/fishpond_results/features")

    # define input paths for features
    dge_path = Path("/data/active/atschan/20220714_apex2_sequencing/fishpond_results/fishpond_results_DMSO_vs_input_residual_lm.csv")
    RBP_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/oRNAment")
    rate_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/Smalec2022/table_s1.xlsx")
    gff_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/fishpond_results/features/gff3/20230724")
    TPM_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/salmon/TPM_salmon_input_DMSO.csv")
    promoter_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/findm/motif_counts")

    # load dge results:
    dge_results = load_dge_results(dge_path)

    # load RBP features
    RBP_features = load_RBP_features(RBP_path)

    # load rate data
    rate_data = load_rate_data(rate_path)

    # load sequence features
    gff_features = load_sequence_features(gff_path)

    # load salmon TPM
    TPM_data = load_TPM(TPM_path)

    # load promoter motifs
    promoter_motifs = load_promoter_motifs(promoter_path)


    # data merging
    df = pd.merge(dge_results[['tx_name', 'log2FC', 'qvalue', 'log10mean', 'SYMBOL','ensembl_transcript_id']], gff_features, on='tx_name')
    df = pd.merge(df, RBP_features, on="ensembl_transcript_id")
    df = pd.merge(df, rate_data, left_on='SYMBOL', right_on='Symbol', how='left')
    df = pd.merge(df, promoter_motifs, left_on='SYMBOL', right_on='SYMBOL')
    df = pd.merge(df, TPM_data[['tx_id', 'TPM']], left_on='ensembl_transcript_id', right_on='tx_id')

    # turn RBP counts into densities
    RBP_cols = RBP_features.drop(['Unnamed: 0', 'ensembl_transcript_id'], axis=1).columns
    df[RBP_cols] = df[RBP_cols].div(df['length_spliced_transcript'], axis=0)
    # turn promoter motif into densities
    promoter_cols = promoter_motifs.columns
    df[promoter_cols] = df[promoter_cols].div(df['length_spliced_transcript'], axis=0)

    #count_cols = [col for col in df if ('_count' in col) & ~('feature_count' in col)]
    #df.drop(count_cols, axis=1, inplace=True)

    # save to csv
    df.to_csv(Path(output_path.joinpath("20230726_transcript_features_merged.csv")))

    # save column names of each feature set separately for later

    column_names = [list(f.columns) for f in [RBP_features, rate_data, gff_features,
                                        TPM_data, promoter_motifs, dge_results]]

    column_names_df = pd.DataFrame({'feature_set': ['RBP_features', 'rate_data', 'gff_features',
                                        'TPM_data', 'promoter_motifs', 'dge_results'],
                                   'column_names': column_names})

    column_names_df.to_csv(output_path.joinpath('20230726_feature_set_columns.csv'))

if __name__ == "__main__":
   main()
