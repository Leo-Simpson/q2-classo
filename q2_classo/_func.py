import numpy as np
from q2_classo.CLasso import *
from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData




def regress(X : FeatureTable[Composition],
            y : MetadataColumn[Numeric],
            #PATH parameters :
            PATH : Bool                         = True,
            PATH_numerical_method : Str         = 'not specified',
            PATH_n_active         : Int         = 0,
            PATH_lambdas          : List[Float] = None,
            PATH_nlam_log         : Int         = 40,
            PATH_lamin_log        : Float       = 1e-2,

            #CV parameters :
            CV : Bool                         = False,
            CV_numerical_method : Str         = 'not specified',
            CV_seed             : Int         = None, # do something here ! for now it can be a bool !
            CV_lambdas          : List[Float] = None, # to do 
            CV_oneSE            : Bool        = True,
            CV_subsets          : Int         = 5,

            #StabSel parameters :
            StabSel : Bool = False,
            StabSel_numerical_method : Str    = 'not specified',
            StabSel_seed             : Int    = None, # do something here ! for now it can be a bool !
            StabSel_lam              : Float  = -1.0, # if negative, then it means 'theoretical'
            StabSel_true_lam         : Bool   = True,
            StabSel_method           : Str    = 'first',
            StabSel_B                : Int    = 50, 
            StabSel_q                : Int    = 10,
            StabSel_percent_nS       : Float  = 0.5,
            StabSel_lamin            : Float  = 1e-2,
            StabSel_threshold        : Float  = 0.7,
            StabSel_threshold_label  : Float  = 0.4, # might unneeded here, but needed for visualisation

            #LAMfixed parameters :
            LAMfixed : Bool = False,
            LAMfixed_numerical_method : Str  = 'not specified',
            LAMfixed_lam              : Float = -1.0, # if negative, then it means 'theoretical'
            LAMfixed_true_lam         : Bool  = True,
            
            #Formulation parameters
            concomitant: Bool      = True,
            huber      : Bool      = False,
            rho        : Float     = 1.345 ,
            rescale    : Bool      = False) -> dict :
    
    X_matrix = np.array(X)
    y_matrix = np.array(y)

    problem = classo_problem(X_matrix, y_matrix, rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho

    if PATH:
        problem.model_selection.PATH= PATH
        param = problem.model_selection.PATHparameters
        param.numerical_method = PATH_numerical_method
        param.n_active         = PATH_n_active
        if PATH_lambdas is None: 
            param.lambdas = np.array([10**(np.log10(PATH_lamin_log) * float(i) / PATH_nlam_log) for i in range(0,PATH_nlam_log) ] )
        else : param.lambdas=  PATH_lambdas

    if CV :
        problem.moodel.selection.LAMfixed = LAMfixed
        param = problem.moodel.selection.CVparameters
        param.numerical_method = CV_numerical_method
        param.seed = CV_numerical_methodCV_seed    
        param.oneSE = CV_numerical_methodCV_oneSE   
        param.Nsubsets = CV_numerical_methodCV_subsets  
        if CV_lambdas is None: param.lambdas =  np.linspace(1., 1e-3, 500)
        else                 : param.lambdas =  CV_lambdas

    if StabSel : 
        problem.moodel.selection.StabSel = StabSel
        param = problem.moodel.selection.StabSelparameters
        param.numerical_method = StabSel_numerical_method
        param.seed = StabSel_seed
        param.true_lam = StabSel_true_lam
        param.method = StabSel_method 
        param.B = StabSel_B
        param.q = StabSel_q
        param.percent_nS = StabSel_percent_nS
        param.lamin = StabSel_lamin
        param.threshold = StabSel_threshold
        param.threshold_label = StabSel_threshold_label
        if (StabSel_lam>0.): param.lam = StabSel_lam
        else               : param.lam = 'theoretical'

    if LAMfixed: 
        problem.moodel.selection.LAMfixed = LAMfixed
        param = problem.moodel.selection.LAMfixedparameters
        param.numerical_method = LAMfixed_numerical_method
        param.true_lam = LAMfixed_true_lam
        if (LAMfixed_lam>0.): param.lam = LAMfixed_lam
        else                : param.lam = 'theoretical'

    problem.solve()

    return problem
    
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





            
