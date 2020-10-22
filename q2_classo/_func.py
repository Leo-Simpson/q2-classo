import numpy as np
import zarr
from classo import *


from qiime2.plugin import (
    SemanticType,
    Plugin,
    Int,
    Float,
    Range,
    Metadata,
    Str,
    Bool,
    Choices,
    MetadataColumn,
    Categorical,
    List,
    Citations,
    TypeMatch,
    Numeric,
)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData
import qiime2
import pandas as pd
import skbio
from ._tree import tree_to_matrix


# TODO : change C type to pandas dataframe


class prediction_data:
    def __init__(self, PATH=False, CV=False, StabSel=False, LAMfixed=False):
        self.PATH = PATH
        self.CV = CV
        self.StabSel = StabSel
        self.LAMfixed = LAMfixed


def generate_data(
    taxa: skbio.TreeNode = None,
    n: int = 100,
    d: int = 80,
    d_nonzero: int = 5,
    classification: bool = False,
) -> (pd.DataFrame, np.ndarray):

    label = np.array(
        ["A" + str(i) for i in range(d // 2)]
        + ["B" + str(i) for i in range(d - d // 2)]
    )
    label_gamma = label[:]
    A = None
    if taxa is not None:
        label2 = np.array([tip.name for tip in taxa.tips()])
        if len(label2) >= d:
            label = label2[:d]
        else:
            label[: len(label2)] = label2
        A, label_gamma = tree_to_matrix(taxa, label)

    (X, C, y), sol = random_data(
        n,
        d,
        d_nonzero,
        0,
        0.5,
        zerosum=True,
        seed=None,
        exp=True,
        A=A,
        classification=classification,
    )

    print(label_gamma[sol != 0])
    if classification:
        y = y == 1.0

    dfx = pd.DataFrame(
        data=X, index=[str(i) for i in range(len(X))], columns=label
    )
    pd.DataFrame(data={"id": range(len(y)), "col": y}).to_csv(
        "randomy.tsv", sep="\t", index=False
    )
    return dfx, C


def transform_features(
    features: pd.DataFrame, transformation: Str = "clr", coef: float = 0.5
) -> pd.DataFrame:
    if transformation == "clr":
        X = features.values
        null_set = X <= 0.0
        X[null_set] = coef
        X = np.log(X)
        X = (X.T - np.mean(X, axis=1)).T

        return pd.DataFrame(
            data=X, index=list(features.index), columns=list(features.columns)
        )

    else:
        raise ValueError(
            "Unknown transformation name, use clr and not %r" % transformation
        )


def add_taxa(
    features: pd.DataFrame, c: np.ndarray = None, taxa: skbio.TreeNode = None
) -> (pd.DataFrame, np.ndarray):

    X = features.values
    label = list(features.columns)
    A, label_new = tree_to_matrix(taxa, label, with_repr=True)
    print(A.shape)
    X_new = X.dot(A)
    if c is None:
        C_new = np.ones((1, len(A))).dot(A)
    else:
        C_new = c.dot(A)
    dfx = pd.DataFrame(
        data=X_new, index=list(features.index), columns=label_new
    )

    return dfx, C_new


def _code_columns(df, column_map, norm=1.0, normalization=False):
    def _code_col(series):
        type_ = column_map[series.index.name].type

        if type_ == "numeric":
            vect = series.to_numpy()
            # if vect has a NaN ==> raise error
            if normalization:
                vect = vect - np.mean(vect)
                vect = vect / np.linalg.norm(vect) * norm
            return vect
        elif type_ == "categorical":
            vect = series.to_numpy()
            # if vect has a NaN ==> raise error
            verfify_binary(vect)
            vect = (
                vect == vect[0]
            )  # set the vector to true if the value is the
            vect = 2 * Y - 1  # transform it to a vector of 1 and -1
            return vect
        else:
            raise ValueError("Unknown type.")

    return df.apply(_code_col, axis=0)


def add_covariates(
    covariates: qiime2.Metadata,
    to_add: str,
    features: pd.DataFrame = None,
    c: np.ndarray = None,
) -> (pd.DataFrame, np.ndarray):

    if features is not None:
        norm = features.apply(np.linalg.norm, axis=0).mean()
        normalization = True
    else:
        norm = 1
        normalization = False

    covariates_df = covariates.to_dataframe()[to_add]
    covariates_df = _code_columns(
        covariate_df, covariates.columns, norm, normalization
    )
    if features is not None:
        # features, covariates_df
        # = features.align(covariates_df, join='inner',axis=0)
        X = pd.concat([features, covariates_df], axis=1, join="inner")

        c_new = np.zeros((len(c), len(c[0]) + len(to_add)))
        c_new[:, : len(c[0])] = c
    else:
        X = covariates_df
        c_new = np.zeros((1, len(to_add)))

    return X, c_new


def regress(
    features: pd.DataFrame,
    y: qiime2.NumericMetadataColumn,
    c: np.ndarray = None,
    do_yshift: bool = False,
    # taxa: skbio.TreeNode = None,
    # PATH parameters :
    path: bool = True,
    path_numerical_method: str = "not specified",
    path_n_active: int = 0,
    path_lambdas: list = None,
    path_nlam_log: int = 40,
    path_lamin_log: float = 1e-2,
    # CV parameters :
    cv: bool = True,
    cv_numerical_method: str = "not specified",
    cv_seed: int = 1,
    cv_lambdas: list = None,  # to do
    cv_one_se: bool = True,
    cv_subsets: int = 5,
    # StabSel parameters :
    stabsel: bool = True,
    stabsel_numerical_method: str = "not specified",
    stabsel_seed: int = None,  # do something here ! for now it can be a bool !
    stabsel_lam: float = -1.0,  # if negative, then it means 'theoretical'
    stabsel_true_lam: bool = True,
    stabsel_method: str = "first",
    stabsel_b: int = 50,
    stabsel_q: int = 10,
    stabsel_percent_ns: float = 0.5,
    stabsel_lamin: float = 1e-2,
    stabsel_threshold: float = 0.7,
    stabsel_threshold_label: float = 0.4,
    # might unneeded here, but needed for visualisation
    # LAMfixed parameters :
    lamfixed: bool = True,
    lamfixed_numerical_method: str = "not specified",
    lamfixed_lam: float = -1.0,  # if negative, then it means 'theoretical'
    lamfixed_true_lam: bool = True,
    # Formulation parameters
    concomitant: bool = True,
    huber: bool = False,
    rho: float = 1.345,
    intercept: bool = True,
) -> classo_problem:

    complete_y = y.to_series()
    complete_y = complete_y[~complete_y.isna()]

    features, pdY = features.align(y.to_series(), join="inner", axis=0)
    missing = pdY.isna()
    training_labels = list(pdY[~missing].index)
    label_missing = list(pdY.index[missing])
    if label_missing:
        print("{} are missing in y ".format(label_missing))
    Y = pdY[~missing].to_numpy()
    X = features.values[~missing, :]

    print(Y.shape, X.shape)

    if do_yshift:
        Y = Y - np.mean(Y)

    problem = classo_problem(X, Y, C=c, label=list(features.columns))
    problem.formulation.huber = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho = rho
    problem.formulation.intercept = intercept

    problem.model_selection.PATH = path
    if path:
        param = problem.model_selection.PATHparameters
        param.numerical_method = path_numerical_method
        param.n_active = path_n_active
        if path_lambdas is None:
            param.lambdas = np.array(
                [
                    10 ** (np.log10(path_lamin_log) * float(i) / path_nlam_log)
                    for i in range(0, path_nlam_log)
                ]
            )
        else:
            param.lambdas = path_lambdas

    problem.model_selection.CV = cv
    if cv:
        param = problem.model_selection.CVparameters
        param.numerical_method = cv_numerical_method
        param.seed = cv_seed
        param.oneSE = cv_one_se
        param.Nsubsets = cv_subsets
        if cv_lambdas is None:
            param.lambdas = np.linspace(1.0, 1e-3, 500)
        else:
            param.lambdas = cv_lambdas

    problem.model_selection.StabSel = stabsel
    if stabsel:
        param = problem.model_selection.StabSelparameters
        param.numerical_method = stabsel_numerical_method
        param.seed = stabsel_seed
        param.true_lam = stabsel_true_lam
        param.method = stabsel_method
        param.B = stabsel_b
        param.q = stabsel_q
        param.percent_nS = stabsel_percent_ns
        param.lamin = stabsel_lamin
        param.threshold = stabsel_threshold
        param.threshold_label = stabsel_threshold_label
        if stabsel_lam > 0.0:
            param.lam = stabsel_lam
        else:
            param.lam = "theoretical"

    problem.model_selection.LAMfixed = lamfixed
    if lamfixed:
        param = problem.model_selection.LAMfixedparameters
        param.numerical_method = lamfixed_numerical_method
        param.true_lam = lamfixed_true_lam
        if lamfixed_lam > 0.0:
            param.lam = lamfixed_lam
        else:
            param.lam = "theoretical"

    problem.solve()

    problem.data.complete_y = complete_y.values
    problem.data.complete_labels = list(complete_y.index)
    problem.data.training_labels = training_labels

    return problem


def classify(
    features: pd.DataFrame,
    y: qiime2.CategoricalMetadataColumn,
    c: np.ndarray = None,
    # taxa: skbio.TreeNode = None,
    # PATH parameters :
    path: bool = True,
    path_numerical_method: str = "not specified",
    path_n_active: int = 0,
    path_lambdas: list = None,
    path_nlam_log: int = 40,
    path_lamin_log: float = 1e-2,
    # CV parameters :
    cv: bool = True,
    cv_numerical_method: str = "not specified",
    cv_seed: int = 1,
    cv_lambdas: list = None,  # to do
    cv_one_se: bool = True,
    cv_subsets: int = 5,
    # StabSel parameters :
    stabsel: bool = True,
    stabsel_numerical_method: str = "not specified",
    stabsel_seed: int = None,  # do something here ! for now it can be a bool !
    stabsel_lam: float = -1.0,  # if negative, then it means 'theoretical'
    stabsel_true_lam: bool = True,
    stabsel_method: str = "first",
    stabsel_b: int = 50,
    stabsel_q: int = 10,
    stabsel_percent_ns: float = 0.5,
    stabsel_lamin: float = 1e-2,
    stabsel_threshold: float = 0.7,
    stabsel_threshold_label: float = 0.4,
    # might unneeded here, but needed for visualisation
    # LAMfixed parameters :
    lamfixed: bool = True,
    lamfixed_numerical_method: str = "not specified",
    lamfixed_lam: float = -1.0,  # if negative, then it means 'theoretical'
    lamfixed_true_lam: bool = True,
    # Formulation parameters
    huber: bool = False,
    rho: float = -1.0,
    intercept: bool = True,
) -> classo_problem:

    features, pdY = features.align(y.to_series(), join="inner", axis=0)
    missing = pdY.isna()
    label_missing = list(pdY.index[missing])
    if label_missing:
        print("{} are missing in y ".format(label_missing))
    Y = pdY[~missing].to_numpy()
    X = features.values[~missing, :]

    verfify_binary(Y)
    Y = Y == Y[0]
    Y = 2 * Y - 1

    problem = classo_problem(X, Y, C=c, label=list(features.columns))
    problem.formulation.classification = True
    problem.formulation.concomitant = False
    problem.formulation.huber = huber
    problem.formulation.rho_classification = rho
    problem.formulation.intercept = intercept

    problem.model_selection.PATH = path
    if path:
        param = problem.model_selection.PATHparameters
        param.numerical_method = path_numerical_method
        param.n_active = path_n_active
        if path_lambdas is None:
            param.lambdas = np.array(
                [
                    10 ** (np.log10(path_lamin_log) * float(i) / path_nlam_log)
                    for i in range(0, path_nlam_log)
                ]
            )
        else:
            param.lambdas = path_lambdas

    problem.model_selection.CV = cv
    if cv:
        param = problem.model_selection.CVparameters
        param.numerical_method = cv_numerical_method
        param.seed = cv_seed
        param.oneSE = cv_one_se
        param.Nsubsets = cv_subsets
        if cv_lambdas is None:
            param.lambdas = np.linspace(1.0, 1e-3, 500)
        else:
            param.lambdas = cv_lambdas

    problem.model_selection.StabSel = stabsel
    if stabsel:
        param = problem.model_selection.StabSelparameters
        param.numerical_method = stabsel_numerical_method
        param.seed = stabsel_seed
        param.true_lam = stabsel_true_lam
        param.method = stabsel_method
        param.B = stabsel_b
        param.q = stabsel_q
        param.percent_nS = stabsel_percent_ns
        param.lamin = stabsel_lamin
        param.threshold = stabsel_threshold
        param.threshold_label = stabsel_threshold_label
        if stabsel_lam > 0.0:
            param.lam = stabsel_lam
        else:
            param.lam = "theoretical"

    problem.model_selection.LAMfixed = lamfixed
    if lamfixed:
        param = problem.model_selection.LAMfixedparameters
        param.numerical_method = lamfixed_numerical_method
        param.true_lam = lamfixed_true_lam
        if lamfixed_lam > 0.0:
            param.lam = lamfixed_lam
        else:
            param.lam = "theoretical"

    problem.solve()

    return problem


def verfify_binary(y):
    lis = []
    for i in y:
        if i not in lis:
            lis.append(i)
    if len(lis) > 2:
        raise ValueError(
            "Metadata column y is supposed to be binary, "
            + "but takes more than 2 different values : "
            + " ; ".join(lis)
        )


def predict(
    features: pd.DataFrame, problem: zarr.hierarchy.Group
) -> prediction_data:
    labels = np.array(problem["data/label"])
    n, d = features.shape[0], len(labels)
    X = np.zeros((n, d))
    for i, label in enumerate(labels):
        if label == "intercept":
            X[:, i] = np.ones(n)
        else:
            X[:, i] = features[label]

    classification = problem["formulation"].attrs["classification"]
    dico_ms = problem["model_selection"].attrs.asdict()
    predictions = prediction_data(
        PATH=problem["model_selection"].attrs["PATH"],
        CV=problem["model_selection"].attrs["CV"],
        StabSel=problem["model_selection"].attrs["StabSel"],
        LAMfixed=problem["model_selection"].attrs["LAMfixed"],
    )

    predictions.sample_labels = list(features.index)

    if predictions.PATH:
        beta_path = np.array(problem["solution/PATH/BETAS"])
        nlam = len(beta_path)
        Y_path = np.zeros((nlam, n))
        for i in range(nlam):
            Y_path[i, :] = X.dot(beta_path[i])
        predictions.YhatPATH = Y_path

    if predictions.CV:
        beta_refit = np.array(problem["solution/CV/refit"])
        Yhat = X.dot(beta_refit)
        predictions.YhatCV = Yhat

    if predictions.StabSel:
        beta_refit = np.array(problem["solution/StabSel/refit"])
        Yhat = X.dot(beta_refit)
        predictions.YhatStabSel = Yhat

    if predictions.LAMfixed:
        beta = np.array(problem["solution/LAMfixed/beta"])
        beta_refit = np.array(problem["solution/LAMfixed/refit"])
        Yhat, Yhatrefit = X.dot(beta), X.dot(beta_refit)
        predictions.YhatLAMfixed = Yhat
        predictions.YhatrefitLAMfixed = Yhatrefit

    return predictions
