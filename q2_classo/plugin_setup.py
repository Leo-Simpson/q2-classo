from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch)


# from q2_classo import regress, classify

#citations = Citations.load('citations.bib', package='q2_classo') 
plugin = Plugin(
name='classo',
                version='0.0.0.dev0',
website='https://github.com/Leo-Simpson/q2-classo',
package='q2-classo',
short_description=('Package for constrained sparse regression and classification'),
description=('This is QIIME 2 plugin that enables sparse and robust linear regression and classification with linear equality constraints on the model parameters.')
)

# remove the - and put _ instead

# plugin.methods.register_function(
#            function=regress,
#            inputs={
#                'X': FeatureTable[Composition],
#                                  },
#            parameters=regress_parameters,
#            outputs=[
#                     # Output about PATH computation
#                     ('PATH_betas',FeatureData[Float]),
#                     ('PATH_sigmas',FeatureData[Float]),
#                     ('PATH_lambdas',FeatureData[Float]),
#                     ('PATH_method',Str),
#                     ('PATH_formulation',Str)
#                     ('PATH_time',Float)
                    
#                     # Output about CV computation
                    
#                     # Output about StabSel computation
                    
#                     # Output about LAMfixed computation
                    
                    
#                     ],
#            input_descriptions=regress_input_descriptions
#            parameter_descriptions=regress_parameter_descriptions,
#            output_descriptions=regress_output_descriptions,
#            name='regress',
#            description=("The function computes the constrainted_sparse_regression vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given"),
#            citations=[citations['Weiss2017']]
#            )





# plugin.methods.register_function(
#                                  function=classify,
#                                  inputs={
#                                  'X': FeatureData[Float],
#                                  'y': FeatureData[Float],
#                                  'C': FeatureData[Float]
#                                  },
#                                  parameters=classify_parameters,
#                                  outputs=[('name',type)],
#                                  input_descriptions=classify_input_descriptions
#                                  parameter_descriptions=classify_parameter_descriptions,
#                                  output_descriptions={
#                                  'rarefied_table': 'The resulting rarefied feature table.'
#                                  },
#                                  name='classify',
#                                  description=("The function computes the constrainted_sparse_linear_classification vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given"),
#                                  citations=[citations['Weiss2017']]
#                                  )









# regress_input_descriptions={
#     'X': 'Matrix representing the data of the problem',
#     'y': 'Vector representing the output of the problem',
#     'C': 'Matrix of constraints to the problem. If it is ‘zero_sum’ then the corresponding attribute will be all_one matrix.'
# }
# regress_parameters={
#      'y': MetadataColumn[Numeric],
#     #Formulation parameters
#     'concomitant': Bool,
#     'huber'      : Bool,
#     'rho'        : Float,
#     'rescale'    : Bool,


#     #PATH parameters :
#     'PATH' : Bool,
#     'PATH_numerical_method' : Str,
#     'PATH_n_active'         : Int, # do something here ! for now it can be a bool !
#     'PATH_lambdas'          : FeatureData[Float],


#     #CV parameters :
#     'CV' : Bool,
#     'CV_numerical_method' : Str,
#     'CV_seed'             : Int, # do something here ! for now it can be a bool !
#     'CV_lambdas'          : FeatureData[Float],
#     'CV_oneSE'            : Bool,
#     'CV_subsets'          : Int,

#     #StabSel parameters :
#     'StabSel' : Bool,
#     'StabSel_numerical_method' : Str,
#     'StabSel_seed'             : Int, # do something here ! for now it can be a bool !
#     'StabSel_lam'              : Float, # can be str as well for now !
#     'StabSel_true_lam'         : Bool,
#     'StabSel_theoretical_lam'  : Float,
#     'StabSel_method'           : Str,
#     'StabSel_B'                : Int,
#     'StabSel_q'                : Int,
#     'StabSel_percent_nS'       : Float,
#     'StabSel_lamin'            : Float,
#     'StabSel_hd'               : Bool,
#     'StabSel_threshold'        : Float,
#     'StabSel_threshold_label'  : Float, # might unneeded here, but needed for visualisation

#     #LAMfixed parameters :
#     'LAMfixed' : Bool,
#     'LAMfixed_numerical_method' : Str,
#     'LAMfixed_lam'              : Float, # can be str as well for now !
#     'LAMfixed_true_lam'         : Bool,
#     'LAMfixed_theoretical_lam'  : Float
# }

# regress_parameter_descriptions={
#     #Formulation parameters
#     'concomitant': 'True if the formulation of the problem should be with an M_estimation of sigma. Default value = True',
#     'huber'      : 'True if the formulation of the problem should be robust Default value = False',
#     'rho'        : 'Value of rho for robust problem. Default value = 1.345',
#     'rescale'    : 'if True, then the function rescale() will be applied to data when solving the problem',


#     #PATH parameters :
#     'PATH' : 'True if path should be computed. Default Value = False',
#     'PATH_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
#     'PATH_n_active'         : 'if it is an integer, then the algo stop computing the path when n_active variables are actives. then the solution does not change from this point. Dafault value : False',
#     'PATH_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.array([10**(-delta * float(i) / nlam) for i in range(0,nlam) ] ) with delta=2. and nlam = 40',

#     #CV parameters :
#     'CV' : 'True if Cross Validation should be computed. Default Value = False',
#     'CV_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
#     'CV_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
#     'CV_lambdas'          : 'list of lambdas for computinf lasso-path for cross validation on lambda. Default value : np.linspace(1., 1e-3, 500)',
#     'CV_oneSE'            : 'if set to True, the selected lambda if computed with method ‘one-standard-error’ Default value : True',
#     'CV_subsets'          : 'number of subset in the cross validation method Dafault value : 5',
#     #StabSel parameters :
#     'StabSel' : 'True if Stability Selection should be computed. Default Value = True',
#     'StabSel_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
#     'StabSel_seed'             : 'Seed for random values, for an equal seed, the result will be the same. If set to False/None: pseudo-random vectors Default value : None',
#     'StabSel_lam'              : '(only used if method = ‘lam’) lam for which the lasso should be computed. Default value : ‘theoretical’ which mean it will be equal to theoretical_lam once it is computed',
#     'StabSel_true_lam'         : '(only used if method = ‘lam’) True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = ‘theoretical’ , then it will takes the value n*theoretical_lam. Default value : True',
#     'StabSel_theoretical_lam'  : '(only used if method = ‘lam’) Theoretical lam. Default value : 0.0 (once it is not computed yet, it is computed thanks to the function theoretical_lam() used in classo_problem.solve())',
#     'StabSel_method'           : '‘first’, ‘lam’ or ‘max’ depending on the type of stability selection we do. Default value : ‘first’',
#     'StabSel_B'                : 'number of subsample considered. Default value : 50',
#     'StabSel_q'                : 'number of selected variable per subsample. Default value : 10',
#     'StabSel_percent_nS'       : 'size of subsample relatively to the total amount of sample Default value : 0.5',
#     'StabSel_lamin'            : 'lamin when computing the lasso-path for method ‘max’ Default value : 1e-2',
#     'StabSel_hd'               : 'if set to True, then the ‘max’ will stop when it reaches n-k actives variables Default value : False',
#     'StabSel_threshold'        : 'threhold for stability selection Default value : 0.7',
#     'StabSel_threshold_label'  : 'threshold to know when the label should be plot on the graph. Default value : 0.4', # might unneeded here, but needed for visualisation

#     #LAMfixed parameters :
#     'LAMfixed' : 'True if solution for a fixed lambda should be computed. Default Value = False',
#     'LAMfixed_numerical_method' : 'name of the numerical method that is used, it can be : ‘Path-Alg’ (path algorithm) , ‘P-PDS’ (Projected primal-dual splitting method) , ‘PF-PDS’ (Projection-free primal-dual splitting method) or ‘DR’ (Douglas-Rachford-type splitting method) Default value : ‘choose’, which means that the function choose_numerical_method() will choose it accordingly to the formulation',
#     'LAMfixed_lam'              : 'lam for which the lasso should be computed. Default value : ‘theoretical’ which mean it will be equal to theoretical_lam once it is computed', # can be str as well for now !
#     'LAMfixed_true_lam'         : 'True if the lambda given is the real lambda, False if it lambda/lambdamax which is between 0 and 1. If True and lam = ‘theoretical’ , then it will takes the value n*theoretical_lam. Default value : True',
#     'LAMfixed_theoretical_lam'  : 'Theoretical lam Default value : 0.0 (once it is not computed yet, it is computed thanks to the function theoretical_lam() used in classo_problem.solve())'

# }

# regress_output_descriptions = {

# }
# }
