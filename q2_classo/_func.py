import numpy as np
from q2_classo.CLasso import *
from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData



Beta_type         = SemanticType('Beta_type') # array of size Npath x d with the solution beta for each lambda on each row
Sigma_type        = SemanticType('Sigma_type') # array of size Npath with the solution sigma for each lambda when the formulation of the problem is R2 or R4
Lambda_type       = SemanticType('Lambda_type') # array of size Npath with the lambdas (real lambdas, not divided by lambda_max) for which the solution is computed
Method_type       = SemanticType('Method_type') # name of the numerical method that has been used. It can be 'Path-Alg', 'P-PDS' , 'PF-PDS' or 'DR'
Formulation_type  = SemanticType('Formulation_type') # can be 'R1' ; 'R2' ; 'R3' ; 'R4' ; 'C1' ; 'C2'
Time_type         = SemanticType('Time_type') # running time of this action

Path_type = SemanticType('Path_type', field_names=['betas','sigmas','lambdas', 'method', 'formulation', 'time'],
                      field_members={
                          'betas' : Beta_type,
                          'sigmas': Sigma_type,
                          'lambdas': Lambda_type,
                          'method' : Method_type,
                          'formulation' : Formulation_type,
                          'time'   : Time_type
                       })

Path = Path_type[Beta_type,Sigma_type,Lambda_type, Method_type, Formulation_type, Time_type]


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
            rescale    : Bool      = False) -> Path_type[Beta_type,Sigma_type,Lambda_type, Method_type, Formulation_type, Time_type]:
    
    
    X_matrix = np.array(X)
    y_matrix = np.array(y)

    problem = classo_problem(X_matrix, y_matrix, rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho

    #PATH
    if PATH:
        problem.model_selection.PATH                            = PATH
        problem.model_selection.PATHparameters.numerical_method = PATH_numerical_method
        problem.model_selection.PATHparameters.n_active         = PATH_n_active
        if PATH_lambdas is None: 
            nlam, lamin = PATH_nlam_log, PATH_lamin_log
            log_lambdas =  np.array([10**(np.log10(lamin) * float(i) / nlam) for i in range(0,nlam) ] )
            problem.model_selection.PATHparameters.lambdas = log_lambdas
        else : problem.model_selection.PATHparameters.lambdas          =  PATH_lambdas

    if CV :
        problem.moodel.selection.LAMfixed = LAMfixed
        if CV_lambdas is None: problem.moodel.selection.CVparameters.lambdas   =  np.linspace(1., 1e-3, 500)
        else                 :  problem.model_selection.CVparameters.lambdas   =  CV_lambdas

    if StabSel : 
        problem.moodel.selection.StabSel = StabSel
        if (StabSel_lam>0.): problem.moodel.selection.StabSelparameters.lam = StabSel_lam
        else               : problem.moodel.selection.StabSelparameters.lam = 'theoretical'

    if LAMfixed: 
        problem.moodel.selection.LAMfixed = LAMfixed
        if (LAMfixed_lam>0.): problem.moodel.selection.LAMfixedparameters.lam = LAMfixed_lam
        else                : problem.moodel.selection.LAMfixedparameters.lam = 'theoretical'


    problem.solve()
    solution_PATH, solution_CV, solution_StabSel, solution_LAM = problem.solution.PATH, problem.solution.CV, problem.solution.StabSel, problem.solution.LAMfixed

    return( solution_PATH.BETAS, solution_PATH.SIGMAS, solution_PATH.LAMBDAS, solution_PATH.method, solution_PATH.formulation, solution_PATH.time  )





            
