import numpy as np
from q2_classo.CLasso import *
from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData
import qiime2



def regress(
            features : np.ndarray,
            y : qiime2.NumericMetadataColumn,
            #PATH parameters :
            path : bool = True,
            path_numerical_method : str         = 'not specified',
            path_n_active         : int         = 0,
            path_lambdas          : list  = None,
            path_nlam_log         : int         = 40,
            path_lamin_log        : float       = 1e-2,

            #CV parameters :
            cv : bool                         = False,
            cv_numerical_method : str         = 'not specified',
            cv_seed             : int         = None, # do something here ! for now it can be a bool !
            cv_lambdas          : list  = None, # to do 
            cv_one_se            : bool        = True,
            cv_subsets          : int         = 5,

            #StabSel parameters :
            stabsel : bool = False,
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
            lamfixed : bool = False,
            lamfixed_numerical_method : str  = 'not specified',
            lamfixed_lam              : float = -1.0, # if negative, then it means 'theoretical'
            lamfixed_true_lam         : bool  = True,
            
            #Formulation parameters
            concomitant: bool      = True,
            huber      : bool      = False,
            rho        : float     = 1.345,
            rescale    : bool      = False) -> classo_problem :


    y = y.to_series().to_numpy()

    if len(y)>len(features): 
        print("More outputs than features ! ")
        y = y[:len(features)]
    elif len(features) > len(y): 
        print("More features than outputs !")
        features = features[:len(y)]


    problem = classo_problem(features, y, rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho

    if path:
        problem.model_selection.PATH= path
        param = problem.model_selection.PATHparameters
        param.numerical_method = path_numerical_method
        param.n_active         = path_n_active
        if path_lambdas is None: 
            param.lambdas = np.array([10**(np.log10(path_lamin_log) * float(i) / path_nlam_log) for i in range(0,path_nlam_log) ] )
        else : param.lambdas=  path_lambdas

    if cv :
        problem.moodel.selection.LAMfixed = LAMfixed
        param = problem.moodel.selection.CVparameters
        param.numerical_method = cv_numerical_method
        param.seed = cv_seed    
        param.oneSE = cv_one_se  
        param.Nsubsets = cv_subsets  
        if cv_lambdas is None: param.lambdas =  np.linspace(1., 1e-3, 500)
        else                 : param.lambdas =  cv_lambdas

    if stabsel : 
        problem.moodel.selection.StabSel = stabsel
        param = problem.moodel.selection.StabSelparameters
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

    if lamfixed: 
        problem.moodel.selection.LAMfixed = lamfixed
        param = problem.moodel.selection.LAMfixedparameters
        param.numerical_method = lamfixed_numerical_method
        param.true_lam = lamfixed_true_lam
        if (lamfixed_lam>0.): param.lam = lamfixed_lam
        else                : param.lam = 'theoretical'

    problem.solve()

    return problem
    
    '''
    solution_PATH, solution_CV, solution_StabSel, solution_LAM = problem.solution.PATH, problem.solution.CV, problem.solution.StabSel, problem.solution.LAMfixed

    output = dict()
    if PATH : output["PATH"] = [solution_PATH.BETAS, solution_PATH.SIGMAS, solution_PATH.LAMBDAS, solution_PATH.method, solution_PATH.formulation, solution_PATH.time]
    else : output["PATH"] = 'not_computed'

    if CV : output["CV"] = [solution_CV.beta]
    else : output["CV"] = 'not_computed'

    if StabSel : output["StabSel"] = [solution_StabSel.beta]
    else : output["StabSel"] = 'not_computed'

    if LAMfixed : output["LAMfixed"] = [solution_LAM.beta]
    else : output["LAMfixed"] = 'not_computed'
    return output
    '''

def generate_data(n : int = 100,
                  d : int = 100,
                  d_nonzero : int = 5
                    ) -> np.ndarray :

    m,d,d_nonzero,k,sigma =100,100,5,1,0.5
    (X,C,y),sol = random_data(n,d,d_nonzero,0,0.5,zerosum=True,seed= 4)
    return X



            
