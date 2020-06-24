import numpy as np
from classo import *


from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData
import qiime2
import pandas as pd
import skbio
from _tree import tree_to_matrix




def generate_data(taxa : skbio.TreeNode = None,
                  n : int = 100,
                  d : int = 80,
                  d_nonzero : int = 5,
                  classification : bool = False
                    ) -> (pd.DataFrame, np.ndarray) :

    
    label = np.array(['A'+str(i) for i in range(d//2)]+['B'+str(i) for i in range(d-d//2)])
    label_gamma = label[:]
    A = None
    if not taxa is None :
        label2 = np.array([tip.name for tip in taxa.tips()] )
        if len(label2 ) >= d : label = label2[:d]
        else : label[:len(label2)] = label2
        A, label_gamma = tree_to_matrix(taxa,label)

    (X,C,y),sol = random_data(n,d,d_nonzero,0,0.5,zerosum=True,seed= None, exp = True, A=A,classification = classification)

    print( label_gamma[ sol != 0 ] )
    if classification : y = y==1.

    dfx = pd.DataFrame(data = X, index = [str(i) for i in range(len(X))] ,columns = label)
    pd.DataFrame(data={'id':range(len(y)),'col':y}).to_csv("randomy.tsv",sep='\t',index=False)
    return dfx, C


def features_clr(features : pd.DataFrame, coef : float = 0.5) -> pd.DataFrame:
    X = features.values
    null_set = (X<=0.)
    X[null_set] = coef
    X = np.log(X)
    X = (X.T - np.mean(X, axis=1)).T

    return pd.DataFrame(data = X, index = list(features.index) ,columns = list(features.columns))


def add_taxa(features : pd.DataFrame,
            c : np.ndarray  = None,
            taxa: skbio.TreeNode = None)-> (pd.DataFrame, np.ndarray) :


    X = features.values
    label = list(features.columns)
    A, label_new = tree_to_matrix(taxa,label, with_repr=True)
    X_new = X.dot(A)
    C_new = c.dot(A)
    dfx = pd.DataFrame(data = X_new, index = list(features.index) ,columns = label_new)
    return dfx, C_new


def add_covariates(features : pd.DataFrame, 
                covariates : qiime2.Metadata, 
                to_add : str,
                c : np.ndarray = None)-> (pd.DataFrame, np.ndarray):

    label = list(features.columns)
    X = features.values
    n,d,k = len(X),len(X.T), len(C)
    


    norm = np.mean( [ np.linalg.norm(X[:,j]) for j in range(d)]  ) # not sure here, maybe set it to 1 instead

    X_new = np.zeros((n,d+len(to_add))  )
    C_new = np.zeros((k,d+len(to_add)))
    X_new[:,:d] = X
    if c is None : C_new[:,:d] = 1.
    else : C_new[:,:d] = c

    covariates = covariates.filter_columns(column_type='numeric')

    for i, name in enumerate(to_add):
        #vect = Y[name].to_numpy() # ???
        try:
            col = covariates.get_column(name)
        except ValueError:
            raise ValueError("Column %r is not numeric or available in the"
                             " metadata" % name)
        vect = col.to_series().to_numpy()
        
        #if category ? 
        #   vect = vect==vect[0] # set the vector to true if the value is the 
        #   vect = 2*vect-1 # transform it to a vector of 1 and -1
        # else : 
        vect  = np.exp(  vect/np.linalg.norm(vect) * norm )   # we take the exp because then in regress or classify, we take the log
        X_new[:,d+i] =  vect
        label.append(name)
    
    dfx = pd.DataFrame(data = X_new, index = [str(i) for i in range(n)] ,columns = label)

    return dfx, C_new


def regress(features : pd.DataFrame,
            y : qiime2.NumericMetadataColumn,
            c : np.ndarray  = None,
            do_yshift : bool = True,
            #taxa: skbio.TreeNode = None,
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
    if do_yshift : Y = Y - np.mean(Y)
    Features = features.values


    problem = classo_problem(Features, Y , C = c, rescale=rescale, label = list(features.columns) )
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
    

def classify(features : pd.DataFrame,
            y : qiime2.CategoricalMetadataColumn,
            c : np.ndarray  = None,
            #taxa: skbio.TreeNode = None,
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
            huber      : bool      = False,
            rho        : float     = -1.,
            rescale    : bool      = False) -> classo_problem :

    y = y.drop_missing_values()
    Y = y.to_series().to_numpy()
    Y = Y==Y[0] # set the vector to true if the value is the 
    Y = 2*Y-1 # transform it to a vector of 1 and -1


    Features = features.values


    problem = classo_problem(Features, Y , C = c, rescale=rescale, label = list(features.columns) )
    problem.formulation.classification = True
    problem.formulation.concomitant = False
    problem.formulation.huber       = huber
    problem.formulation.rho_classification = rho

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
    