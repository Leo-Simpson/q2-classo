import numpy as np
from CLasso import *

def regress(X : FeatureTable[Composition],
            y: MetadataColumn[Numeric],
            #PATH parameters :
            PATH : Bool                                = True,
            PATH_numerical_method : Str                = 'not specified',
            PATH_n_active         : Int                = 0,
            PATH_lambdas          : FeatureData[Float] = None,
            PATH_nlam_log         : Int                = 40,
            PATH_lamin_log        : Float              = 1e-2,

            # #CV parameters :
            # CV : Bool,
            # CV_numerical_method : Str,
            # CV_seed             : Int, # do something here ! for now it can be a bool !
            # CV_lambdas          : FeatureData[Float],
            # CV_oneSE            : Bool,
            # CV_subsets          : Int,

            # #StabSel parameters :
            # StabSel : Bool,
            # StabSel_numerical_method : Str,
            # StabSel_seed             : Int, # do something here ! for now it can be a bool !
            # StabSel_lam              : Float, # can be str as well for now !
            # StabSel_true_lam         : Bool,
            # StabSel_theoretical_lam  : Float,
            # StabSel_method           : Str,
            # StabSel_B                : Int,
            # StabSel_q                : Int,
            # StabSel_percent_nS       : Float,
            # StabSel_lamin            : Float,
            # StabSel_hd               : Bool,
            # StabSel_threshold        : Float,
            # StabSel_threshold_label  : Float, # might unneeded here, but needed for visualisation

            # #LAMfixed parameters :
            # LAMfixed : Bool,
            # LAMfixed_numerical_method : Str,
            # LAMfixed_lam              : Float, # can be str as well for now !
            # LAMfixed_true_lam         : Bool,
            # LAMfixed_theoretical_lam  : Float,
            
            #Formulation parameters
            concomitant: Bool  = True,
            huber      : Bool  = False,
            rho        : Float = 1.345 ,
            rescale    : Bool  = False) -> (
                     FeatureData[Float], FeatureData[Float], FeatureData[Float], Str, Str, Float):
    
    
    problem = classo_problem(np.array(X),np.array(y), rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho
    #PATH
    problem.model_selection.PATH                            = PATH
    problem.model_selection.PATHparameters.numerical_method = PATH_numerical_method
    problem.model_selection.PATHparameters.n_active         = PATH_n_active
    if PATH_lambdas is None: 
        nlam, lamin = PATH_nlam_log, PATH_lamin_log
        log_lambdas =  np.array([10**(np.log10(lamin) * float(i) / nlam) for i in range(0,nlam) ] )
        problem.model_selection.PATHparameters.lambdas = log_lambdas
    else : problem.model_selection.PATHparameters.lambdas          =  PATH_lambdas

    #CV
    #StabSel
    #LAM


    problem.solve()
    solution_PATH, solution_CV, solution_StabSel, solution_LAM = problem.solution.PATH, problem.solution.CV, problem.solution.StabSel, problem.solution.LAMfixed

    return( solution_PATH.BETAS, solution_PATH.SIGMAS, solution_PATH.LAMBDAS, solution_PATH.method, solution_PATH.formulation, solution_PATH.time  )





            
