from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch)


from q2_classo import regress, classify

citations = Citations.load('citations.bib', package='q2_feature_table') # ???
plugin = Plugin(
name='classo',
                version=q2_feature_table.__version__, # ????
website='https://github.com/Leo-Simpson/q2-classo',
package='q2-classo',
short_description=('Package for constrained sparse regression and classification'),
description=('This is QIIME 2 plugin that enables sparse and robust linear regression and classification with linear equality constraints on the model parameters.')
)

plugin.methods.register_function(
           function=regress,
           inputs={
                                 'X': FeatureData[Float],
                                 'y': FeatureData[Float],
                                 'C': FeatureData[Float]
                                 },
           parameters=regress_parameters,
           outputs=[
                    # Output about PATH computation
                    ('PATH-betas',FeatureData[Float]),
                    ('PATH-sigmas',FeatureData[Float]),
                    ('PATH-lambdas',FeatureData[Float]),
                    ('PATH-method',Str),
                    ('PATH-formulation',Str)
                    ('PATH-time',Float)
                    
                    # Output about CV computation
                    
                    # Output about StabSel computation
                    
                    # Output about LAMfixed computation
                    
                    
                    ],
           input_descriptions=regress_input_descriptions
           parameter_descriptions=regress_parameter_descriptions,
           output_descriptions=regress_output_descriptions,
           name='regress',
           description=("The function computes the constrainted-sparse-regression vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given"),
           citations=[citations['Weiss2017']]
           )





plugin.methods.register_function(
                                 function=classify,
                                 inputs={
                                 'X': FeatureData[Float],
                                 'y': FeatureData[Float],
                                 'C': FeatureData[Float]
                                 },
                                 parameters=classify_parameters,
                                 outputs=[('name',type)],
                                 input_descriptions=classify_input_descriptions
                                 parameter_descriptions=classify_parameter_descriptions,
                                 output_descriptions={
                                 'rarefied_table': 'The resulting rarefied feature table.'
                                 },
                                 name='classify',
                                 description=("The function computes the constrainted-sparse-linear-classification vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given"),
                                 citations=[citations['Weiss2017']]
                                 )









regress_input_descriptions={
    'X': 'Matrix representing the data of the problem',
    'y': 'Vector representing the output of the problem',
    'C': 'Matrix of constraints to the problem. If it is ‘zero-sum’ then the corresponding attribute will be all-one matrix.'
}
regress_parameters={
    #Formulation parameters
    'concomitant': Bool,
    'huber'      : Bool,
    'rho'        : Float,
    'rescale'    : Bool,


    #PATH parameters :
    'PATH' : Bool,
    'PATH-numerical_method' : Str,
    'PATH-n_active'         : Int, # do something here ! for now it can be a bool !
    'PATH-lambdas'          : FeatureData[Float],


    #CV parameters :
    'CV' : Bool,
    'CV-numerical_method' : Str,
    'CV-seed'             : Int, # do something here ! for now it can be a bool !
    'CV-lambdas'          : FeatureData[Float],
    'CV-oneSE'            : Bool,
    'CV-subsets'          : Int,

    #StabSel parameters :
    'StabSel' : Bool,
    'StabSel-numerical_method' : Str,
    'StabSel-seed'             : Int, # do something here ! for now it can be a bool !
    'StabSel-lam'              : Float, # can be str as well for now !
    'StabSel-true_lam'         : Bool,
    'StabSel-theoretical_lam'  : Float,
    'StabSel-method'           : Str,
    'StabSel-B'                : Int,
    'StabSel-q'                : Int,
    'StabSel-percent_nS'       : Float,
    'StabSel-lamin'            : Float,
    'StabSel-hd'               : Bool,
    'StabSel-threshold'        : Float,
    'StabSel-threshold-label'  : Float, # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'LAMfixed' : Bool,
    'LAMfixed-numerical_method' : Str,
    'LAMfixed-lam'              : Float, # can be str as well for now !
    'LAMfixed-true_lam'         : Bool,
    'LAMfixed-theoretical_lam'  : Float
}

regress_parameter_descriptions={
    #Formulation parameters
    'concomitant': 'True if the formulation of the problem should be with an M-estimation of sigma. Default value = True',
    'huber'      : 'True if the formulation of the problem should be robust Default value = False',
    'rho'        : 'Value of rho for robust problem. Default value = 1.345',
    'rescale'    : 'if True, then the function rescale() will be applied to data when solving the problem',


    #PATH parameters :
    'PATH' : 'True if path should be computed. Default Value = False',
    'PATH-numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'PATH-n_active'         : 'if it is an integer, then the algo stop computing the path when n_active variables are actives. then the solution does not change from this point. Dafault value : False',
    'PATH-lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.array([10**(-delta * float(i) / nlam) for i in range(0,nlam) ] ) with delta=2. and nlam = 40',

    #CV parameters :
    'CV' : 'True if Cross Validation should be computed. Default Value = False',
    'CV-numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'CV-seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'CV-lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.linspace(1., 1e-3, 500)',
    'CV-oneSE'            : 'if set to True, the selected lambda if computed with method ‘one-standard-error’ Default value : True',
    'CV-subsets'          : 'number of subset in the cross validation method Dafault value : 5',
    #StabSel parameters :
    'StabSel' : 'True if Stability Selection should be computed. Default Value = True',
    'StabSel-numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'StabSel-seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
    'StabSel-lam'              : '(only used if method = ‘lam’) lam for which the lasso should be computed. Default value : ‘theoretical’ which mean it will be equal to theoretical_lam once it is computed',
    'StabSel-true_lam'         : '(only used if method = ‘lam’) True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = ‘theoretical’ , then it will takes the value n*theoretical_lam. Default value : True',
    'StabSel-theoretical_lam'  : '(only used if method = ‘lam’) Theoretical lam. Default value : 0.0 (once it is not computed yet, it is computed thanks to the function theoretical_lam() used in classo_problem.solve())',
    'StabSel-method'           : '‘first’, ‘lam’ or ‘max’ depending on the type of stability selection we do. Default value : ‘first’',
    'StabSel-B'                : 'number of subsample considered. Default value : 50',
    'StabSel-q'                : 'number of selected variable per subsample. Default value : 10',
    'StabSel-percent_nS'       : 'size of subsample relatively to the total amount of sample Default value : 0.5',
    'StabSel-lamin'            : 'lamin when computing the lasso-path for method ‘max’ Default value : 1e-2',
    'StabSel-hd'               : 'if set to True, then the ‘max’ will stop when it reaches n-k actives variables Default value : False',
    'StabSel-threshold'        : 'threhold for stability selection Default value : 0.7',
    'StabSel-threshold-label'  : 'threshold to know when the label should be plot on the graph. Default value : 0.4', # might unneeded here, but needed for visualisation

    #LAMfixed parameters :
    'LAMfixed' : 'True if solution for a fixed lambda should be computed. Default Value = False',
    'LAMfixed-numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
    'LAMfixed-lam'              : 'lam for which the lasso should be computed. Default value : ‘theoretical’ which mean it will be equal to theoretical_lam once it is computed', # can be str as well for now !
    'LAMfixed-true_lam'         : 'True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = ‘theoretical’ , then it will takes the value n*theoretical_lam. Default value : True',
    'LAMfixed-theoretical_lam'  : 'Theoretical lam Default value : 0.0 (once it is not computed yet, it is computed thanks to the function theoretical_lam() used in classo_problem.solve())'

}

regress_output_descriptions = {

}
}
