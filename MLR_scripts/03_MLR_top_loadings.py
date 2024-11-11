from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

from sklearn.model_selection import cross_validate
from sklearn.model_selection import RandomizedSearchCV
from sklearn.linear_model import HuberRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import IsolationForest
from yellowbrick.regressor import CooksDistance
from sklearn.preprocessing import FunctionTransformer
from sklearn.pipeline import Pipeline
from sklearn import set_config

set_config(transform_output="pandas")

def get_top_loadings(pca, column_names, pc):

    loadings = pd.DataFrame(np.abs(pca.components_))
    loadings.columns = column_names
    top_loadings = pd.DataFrame(loadings.iloc[pc, :].sort_values(ascending=False)).iloc[0:10, :].reset_index().rename(columns={'index': 'feature_name', pc: 'loading'})
    top_loadings['loading_rank'] = np.arange(0, 10)
    return top_loadings


def make_regression(df):

    X_train, y_train, X_test, y_test, pca, feature_names = prepare_data(df)

    regr = LinearRegression()
    regr.fit(X_train, y_train)

    # get the top_pcs and their loadings
    top_pcs = np.argsort(np.abs(regr.coef_))[::-1][0:10]

    loadings = []

    for i, pc in enumerate(top_pcs):
        current_top_loadings = get_top_loadings(pca, feature_names, pc)
        current_top_loadings['pc'] = pc
        current_top_loadings['coef'] = regr.coef_[pc]
        current_top_loadings['model_test_score'] = regr.score(X_test, y_test)
        current_top_loadings['pc_coef_rank'] = i
        loadings.append(current_top_loadings)

    top_loadings = pd.concat(loadings)

    return top_loadings

def train_test_split_custom(data, X_pca, train_size):
    train_names = np.random.choice(data.gene_name.unique(), int(len(
        data.gene_name.unique()) * train_size), replace=False)

    train_set_ids = [c for c in
                     data.loc[data.gene_name.isin(train_names)].index]
    test_set_ids = [c for c in
                    data.loc[~data.gene_name.isin(train_names)].index]

    X_pca_df = pd.DataFrame(X_pca, index=data.index)

    X_train = X_pca_df.loc[X_pca_df.index.isin(train_set_ids)]
    X_test = X_pca_df.loc[X_pca_df.index.isin(test_set_ids)]

    y_train = data.loc[data.index.isin(train_set_ids)]['log2FC']
    y_test = data.loc[data.index.isin(test_set_ids)]['log2FC']

    return (X_train, y_train, X_test, y_test)


def prepare_data(df):

    data = df.copy()
    data = data.loc[data.gene_type == 'protein_coding']
    # linear Regression and PCA can't deal with NA, so we need to impute them. Here by taking the mean of the column
    features_for_training = [col for col in data.columns if
                             col not in ['hit_type', 'hit_type_training',
                                         'hit_type_categorical', 'log2FC',
                                         'qvalue',
                                         'Unnamed: 0', 'tx_name', 'gene_id',
                                         'SYMBOL', 'tx_type', 'gene_type',
                                         'gene_id',
                                         'level_1', 'gene_name',
                                         'ensembl_transcript_id',
                                         'Unnamed: 0_x', 'Unnamed: 0_y',
                                         'log2FoldChange', 'fpkm',
                                         'replicate', 'cluster', 'padj',
                                         'descriptor', 'pvalue', 'stat',
                                         'lfcSE', 'tx_id', 'Gene name',
                                         'log10mean']]

    X = data[features_for_training]

    X.replace([np.inf, -np.inf, 0], np.nan, inplace=True)

    # remove columns where > 75% of values are 0
    X.dropna(thresh=X.shape[0]*0.25, axis=1, inplace=True)

    # remove outliers
    y = data['log2FC']

    X_train, y_train, X_test, y_test = train_test_split_custom(data, X, 0.85)

    pipe = Pipeline([('transformer', FunctionTransformer(np.log10)),
                     ('imputer', SimpleImputer()),
                     ('scaler', StandardScaler()),
                     ('pca', PCA(n_components=0.99))])

    X_train = pipe.fit_transform(X_train)
    X_test = pipe.transform(X_test)

    return X_train, y_train, X_test, y_test, pipe['pca'], pipe.feature_names_in_



def main():

    # define output path
    output_path = Path(
        r"/data/active/atschan/20220714_apex2_sequencing/MLR_results/20230726")

    # define path to features
    feature_path = Path(
        r"/data/active/atschan/20220714_apex2_sequencing/fishpond_results/features/20230726_transcript_features_merged.csv")

    # load features
    df = pd.read_csv(feature_path)

    top_loadings_list = []
    for i in range(0, 100):
        top_loadings = make_regression(df)
        top_loadings['iteration'] = i
        top_loadings_list.append(top_loadings)

    top_loadings_df = pd.concat(top_loadings_list)
    top_loadings_df.to_csv(Path(output_path.joinpath("top_loadings_df.csv")))


if __name__ == "__main__":
   main()
