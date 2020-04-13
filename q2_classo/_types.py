from qiime2.plugin import SemanticType
from CLasso import *
from q2_classo.plugin_setup import plugin
from ._formats import CLASSOProblemDirectoryFormat

CLASSOProblem = SemanticType("CLASSOProblem")

ConstraintMatrix = SemanticType("ConstraintMatrix")

plugin.register_semantic_types(CLASSOProblem)

plugin.register_semantic_type_to_format(CLASSOProblem, 
                                        artifact_format=CLASSOProblemDirectoryFormat)





'''
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
'''