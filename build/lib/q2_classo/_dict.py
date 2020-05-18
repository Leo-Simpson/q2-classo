from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition


regress_parameters={
    'y': MetadataColumn[Numeric],
    'do_clr' : Bool,
    #Formulation parameters
    'concomitant': Bool,
    'huber'      : Bool,
    'rho'        : Float,
    'rescale'    : Bool,


    #PATH parameters :
    'path' : Bool,
    'path_numerical_method' : Str,
    'path_n_active'         : Int, # do something here ! for now it can be a bool !
    'path_lambdas'          : List[Float],
    'path_nlam_log'         : Int,
    'path_lamin_log'        : Float,


    #CV parameters :
    'cv' : Bool,
    'cv_numerical_method' : Str,
    'cv_seed'             : Int, # do something here ! for now it can be a bool !
    'cv_lambdas'          : List[Float],
    'cv_one_se'            : Bool,
    'cv_subsets'          : Int,

    #StabSel parameters :
    'stabsel' : Bool,
    'stabsel_numerical_method' : Str,
    'stabsel_seed'             : Int, # do something here ! for now it can be a bool !
    'stabsel_lam'              : Float, # can be str as well for now !
    'stabsel_true_lam'         : Bool,
    'stabsel_method'           : Str,
    'stabsel_b'                : Int,
    'stabsel_q'                : Int,
    'stabsel_percent_ns'       : Float,
    'stabsel_lamin'            : Float,
    'stabsel_threshold'        : Float,
    'stabsel_threshold_label'  : Float, # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'lamfixed' : Bool,
    'lamfixed_numerical_method' : Str,
    'lamfixed_lam'              : Float, # can be str as well for now !
    'lamfixed_true_lam'         : Bool
}
regress_parameter_descriptions={
    'y': 'Vector representing the output of the problem',
    'do_clr' : 'if set to true, then features will be centered-log-ration transformed and y will be centered',
    #Formulation parameters
    'concomitant': 'True if the formulation of the problem should be with an M_estimation of sigma. Default value = True',
    'huber'      : 'True if the formulation of the problem should be robust Default value = False',
    'rho'        : 'Value of rho for robust problem. Default value = 1.345',
    'rescale'    : 'if True, then the function rescale() will be applied to data when solving the problem',


    #PATH parameters :
    'path' : 'True if path should be computed. Default Value = False',
    'path_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'path_n_active'         : 'if it is an integer, then the algo stop computing the path when n_active variables are actives. then the solution does not change from this point. Dafault value : False',
    'path_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.array([10**(-delta * float(i) / nlam) for i in range(0,nlam) ] ) with delta=2. and nlam = 40',
    'path_nlam_log'         : ' number of lambdas required, if the list of lambdas is not specified, in order to use a log-ratio list of lambdas',
    'path_lamin_log'        : 'minimum of lambdas required, if the list of lambdas is not specified, in order to use a log-ratio list of lambdas',

    #CV parameters :
    'cv' : 'True if Cross Validation should be computed. Default Value = False',
    'cv_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'cv_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'cv_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.linspace(1., 1e-3, 500)',
    'cv_one_se'            : 'if set to True, the selected lambda if computed with method ‘one-standard-error’ Default value : True',
    'cv_subsets'          : 'number of subset in the cross validation method Dafault value : 5',
    #StabSel parameters :
    'stabsel' : 'True if Stability Selection should be computed. Default Value = True',
    'stabsel_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'stabsel_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'stabsel_lam'              : '(only used if method = ‘lam’) lam for which the lasso should be computed. Default value : -1. which means it will be equal to theoretical_lam once it is computed',
    'stabsel_true_lam'         : '(only used if method = ‘lam’) True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = -1. , then it will takes the value n*theoretical_lam. Default value : True',
    'stabsel_method'           : '‘first’, ‘lam’ or ‘max’ depending on the type of stability selection we do. Default value : ‘first’',
    'stabsel_b'                : 'number of subsample considered. Default value : 50',
    'stabsel_q'                : 'number of selected variable per subsample. Default value : 10',
    'stabsel_percent_ns'       : 'size of subsample relatively to the total amount of sample Default value : 0.5',
    'stabsel_lamin'            : 'lamin when computing the lasso-path for method ‘max’ Default value : 1e-2',
    'stabsel_threshold'        : 'threhold for stability selection Default value : 0.7',
    'stabsel_threshold_label'  : 'threshold to know when the label should be plot on the graph. Default value : 0.4', # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'lamfixed' : 'True if solution for a fixed lambda should be computed. Default Value = False',
    'lamfixed_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'lamfixed_lam'              : 'lam for which the lasso should be computed. Default value : -1 which mean it will be equal to theoretical_lam once it is computed', # can be str as well for now !
    'lamfixed_true_lam'         : 'True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = -1. , then it will takes the value n*theoretical_lam. Default value : True'
}



classify_parameters={
    'y': MetadataColumn[Categorical],
    'do_clr' : Bool,
    #Formulation parameters
    'huber'      : Bool,
    'rho'        : Float,
    'rescale'    : Bool,


    #PATH parameters :
    'path' : Bool,
    'path_numerical_method' : Str,
    'path_n_active'         : Int, # do something here ! for now it can be a bool !
    'path_lambdas'          : List[Float],
    'path_nlam_log'         : Int,
    'path_lamin_log'        : Float,


    #CV parameters :
    'cv' : Bool,
    'cv_numerical_method' : Str,
    'cv_seed'             : Int, # do something here ! for now it can be a bool !
    'cv_lambdas'          : List[Float],
    'cv_one_se'            : Bool,
    'cv_subsets'          : Int,

    #StabSel parameters :
    'stabsel' : Bool,
    'stabsel_numerical_method' : Str,
    'stabsel_seed'             : Int, # do something here ! for now it can be a bool !
    'stabsel_lam'              : Float, # can be str as well for now !
    'stabsel_true_lam'         : Bool,
    'stabsel_method'           : Str,
    'stabsel_b'                : Int,
    'stabsel_q'                : Int,
    'stabsel_percent_ns'       : Float,
    'stabsel_lamin'            : Float,
    'stabsel_threshold'        : Float,
    'stabsel_threshold_label'  : Float, # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'lamfixed' : Bool,
    'lamfixed_numerical_method' : Str,
    'lamfixed_lam'              : Float, # can be str as well for now !
    'lamfixed_true_lam'         : Bool
}
classify_parameter_descriptions={
    'y': 'Vector representing the output of the problem',
    'do_clr' : 'if set to true, then features will be centered-log-ration transformed and y will be centered',
    #Formulation parameters
    'huber'      : 'True if the formulation of the problem should be robust Default value = False',
    'rho'        : 'Value of rho for robust problem. Default value = 1.345',
    'rescale'    : 'if True, then the function rescale() will be applied to data when solving the problem',


    #PATH parameters :
    'path' : 'True if path should be computed. Default Value = False',
    'path_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'path_n_active'         : 'if it is an integer, then the algo stop computing the path when n_active variables are actives. then the solution does not change from this point. Dafault value : False',
    'path_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.array([10**(-delta * float(i) / nlam) for i in range(0,nlam) ] ) with delta=2. and nlam = 40',
    'path_nlam_log'         : ' number of lambdas required, if the list of lambdas is not specified, in order to use a log-ratio list of lambdas',
    'path_lamin_log'        : 'minimum of lambdas required, if the list of lambdas is not specified, in order to use a log-ratio list of lambdas',

    #CV parameters :
    'cv' : 'True if Cross Validation should be computed. Default Value = False',
    'cv_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'cv_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'cv_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.linspace(1., 1e-3, 500)',
    'cv_one_se'            : 'if set to True, the selected lambda if computed with method ‘one-standard-error’ Default value : True',
    'cv_subsets'          : 'number of subset in the cross validation method Dafault value : 5',
    #StabSel parameters :
    'stabsel' : 'True if Stability Selection should be computed. Default Value = True',
    'stabsel_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'stabsel_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'stabsel_lam'              : '(only used if method = ‘lam’) lam for which the lasso should be computed. Default value : -1. which mean it will be equal to theoretical_lam once it is computed',
    'stabsel_true_lam'         : '(only used if method = ‘lam’) True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = -1. , then it will takes the value n*theoretical_lam. Default value : True',
    'stabsel_method'           : '‘first’, ‘lam’ or ‘max’ depending on the type of stability selection we do. Default value : ‘first’',
    'stabsel_b'                : 'number of subsample considered. Default value : 50',
    'stabsel_q'                : 'number of selected variable per subsample. Default value : 10',
    'stabsel_percent_ns'       : 'size of subsample relatively to the total amount of sample Default value : 0.5',
    'stabsel_lamin'            : 'lamin when computing the lasso-path for method ‘max’ Default value : 1e-2',
    'stabsel_threshold'        : 'threhold for stability selection Default value : 0.7',
    'stabsel_threshold_label'  : 'threshold to know when the label should be plot on the graph. Default value : 0.4', # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'lamfixed' : 'True if solution for a fixed lambda should be computed. Default Value = False',
    'lamfixed_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'lamfixed_lam'              : 'lam for which the lasso should be computed. Default value : -1. which mean it will be equal to theoretical_lam once it is computed', # can be str as well for now !
    'lamfixed_true_lam'         : 'True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = -1. , then it will takes the value n*theoretical_lam. Default value : True'
}




