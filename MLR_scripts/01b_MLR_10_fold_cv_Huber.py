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


def make_regression(df):

    test_scores = []
    X_fold_list, y_fold_list = prepare_data(df)

    for i in range(0, 10):
        # define train and test set
        X_test = X_fold_list[i]
        y_test = y_fold_list[i]

        X_train = X_fold_list.copy()
        X_train.pop(i)
        X_train = pd.concat(X_train)
        y_train = y_fold_list.copy()
        y_train.pop(i)
        y_train = pd.concat(y_train)

        regr = HuberRegressor(max_iter=500)
        regr.fit(X_train, y_train)

        test_scores.append(regr.score(X_test, y_test))

    cv_result = pd.DataFrame({'iteration': np.arange(0, 10),
                              'test_score': test_scores})

    return cv_result

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

def Kfold_custom(data, X, n_splits=10):
    X_df = pd.DataFrame(X, index=data.index)
    # to do: shuffle names
    gene_names = data.gene_name.unique()
    np.random.shuffle(gene_names)
    names_per_fold = np.array_split(gene_names, n_splits)

    X_fold_list = []
    y_fold_list = []

    for fold in names_per_fold:
        ids = [c for c in
                         data.loc[data.gene_name.isin(fold)].index]
        X_fold_list.append(X_df.loc[X_df.index.isin(ids)])
        y_fold_list.append(data.loc[data.index.isin(ids)]['log2FC'])

    return X_fold_list, y_fold_list



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
    X.dropna(thresh=X.shape[0] * 0.25, axis=1, inplace=True)

    # remove outliers
    y = data['log2FC']

    # save X for later

    pipe = Pipeline([('transformer', FunctionTransformer(np.log10)),
                     ('imputer', SimpleImputer()),
                     ('scaler', StandardScaler())])

    X_transformed_imputed_scaled = pipe.fit_transform(X)
    X_transformed_imputed_scaled = pd.DataFrame(X_transformed_imputed_scaled, columns=pipe['imputer'].get_feature_names_out())

    X_transformed_imputed_scaled.to_csv(Path(r"/data/active/atschan/20220714_apex2_sequencing/MLR_results/20230726/feature_input_MLR_preprocessed.csv"))

    # continue with MLR preprocessing
    #X_train, y_train, X_test, y_test = train_test_split_custom(data, X, 0.85)

    pipe = Pipeline([('transformer', FunctionTransformer(np.log10)),
                     ('imputer', SimpleImputer()),
                     ('scaler', StandardScaler()),
                     ('pca', PCA(n_components=0.99))])

    #X_train = pipe.fit_transform(X_train)
    #X_test = pipe.transform(X_test)

    X_processed = pipe.fit_transform(X)

    X_fold_list, y_fold_list = Kfold_custom(data, X_processed, n_splits=10)

    return X_fold_list, y_fold_list



def main():

    # define output path
    output_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/MLR_results/20230726")

    # define path to features
    feature_path = Path(r"/data/active/atschan/20220714_apex2_sequencing/fishpond_results/features/20230726_transcript_features_merged.csv")

    # load features
    df = pd.read_csv(feature_path)
    print('loaded df')

    # run MLR
    cv_result = make_regression(df)

    cv_result.to_csv(output_path.joinpath('20230726_cv_result_Huber.csv'))


if __name__ == "__main__":
   main()
