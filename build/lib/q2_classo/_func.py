import numpy as np
from q2_classo.CLasso import *
from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData
import qiime2
import pandas as pd



def regress(
            features : np.ndarray,
            y : qiime2.NumericMetadataColumn,
            c : np.ndarray  = None,
            do_clr : bool = True,
            #phylogenetic_tree : np.ndarray = None,
            #PATH parameters :
            path : bool = True,
            path_numerical_method : str         = 'not specified',
            path_n_active         : int         = 0,
            path_lambdas          : list  = None,
            path_nlam_log         : int         = 40,
            path_lamin_log        : float       = 1e-2,

            #CV parameters :
            cv : bool                         = True,
            cv_numerical_method : str         = 'not specified',
            cv_seed             : int         = 1,
            cv_lambdas          : list  = None, # to do 
            cv_one_se            : bool        = True,
            cv_subsets          : int         = 5,

            #StabSel parameters :
            stabsel : bool = True,
            stabsel_numerical_method : str    = 'not specified',
            stabsel_seed             : int    = None, # do something here ! for now it can be a bool !
            stabsel_lam              : float  = -1.0, # if negative, then it means 'theoretical'
            stabsel_true_lam         : bool   = True,
            stabsel_method           : str    = 'first',
            stabsel_b                : int    = 50, 
            stabsel_q                : int    = 10,
            stabsel_percent_ns       : float  = 0.5,
            stabsel_lamin            : float  = 1e-2,
            stabsel_threshold        : float  = 0.7,
            stabsel_threshold_label  : float  = 0.4, # might unneeded here, but needed for visualisation

            #LAMfixed parameters :
            lamfixed : bool = True,
            lamfixed_numerical_method : str  = 'not specified',
            lamfixed_lam              : float = -1.0, # if negative, then it means 'theoretical'
            lamfixed_true_lam         : bool  = True,
            
            #Formulation parameters
            concomitant: bool      = True,
            huber      : bool      = False,
            rho        : float     = 1.345,
            rescale    : bool      = False) -> classo_problem :


    Y = y.to_series().to_numpy()
    if do_clr : 
        Features = clr(features.T).T
        Y = Y - np.mean(Y)
    else : 
        Features = features


    problem = classo_problem(Features, Y , C = c, rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho

    problem.model_selection.PATH= path
    if path:
        param = problem.model_selection.PATHparameters
        param.numerical_method = path_numerical_method
        param.n_active         = path_n_active
        if path_lambdas is None: 
            param.lambdas = np.array([10**(np.log10(path_lamin_log) * float(i) / path_nlam_log) for i in range(0,path_nlam_log) ] )
        else : param.lambdas=  path_lambdas

    problem.model_selection.CV = cv
    if cv :
        param = problem.model_selection.CVparameters
        param.numerical_method = cv_numerical_method
        param.seed = cv_seed    
        param.oneSE = cv_one_se  
        param.Nsubsets = cv_subsets  
        if cv_lambdas is None: param.lambdas =  np.linspace(1., 1e-3, 500)
        else                 : param.lambdas =  cv_lambdas

    problem.model_selection.StabSel = stabsel
    if stabsel : 
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
        if (stabsel_lam>0.): param.lam = stabsel_lam
        else               : param.lam = 'theoretical'

    problem.model_selection.LAMfixed = lamfixed
    if lamfixed: 
        param = problem.model_selection.LAMfixedparameters
        param.numerical_method = lamfixed_numerical_method
        param.true_lam = lamfixed_true_lam
        if (lamfixed_lam>0.): param.lam = lamfixed_lam
        else                : param.lam = 'theoretical'

    problem.solve()

    return problem
    

def generate_data(n : int = 100,
                  d : int = 150,
                  d_nonzero : int = 5
                    ) -> (np.ndarray, np.ndarray) :

    (X,C,y),sol = random_data(n,d,d_nonzero,0,0.5,zerosum=True,seed= 4, exp = True)
    pd.DataFrame(data={'id':range(len(y)),'col':y}).to_csv("randomy.tsv",sep='\t',index=False)
    return X, C

