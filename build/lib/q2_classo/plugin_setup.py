from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData

from . import  *



#citations = Citations.load('citations.bib', package='q2_classo') 
plugin = Plugin(
name='classo',
                version='0.0.0.dev0',
website='https://github.com/Leo-Simpson/q2-classo',
package='q2-classo',
short_description=('Package for constrained sparse regression and classification'),
description=('This is QIIME 2 plugin that enables sparse and robust linear regression and classification with linear equality constraints on the model parameters.')
)


plugin.methods.register_function(
           function=regress,
           inputs={'X': FeatureTable[Composition]},
           parameters=regress_parameters,
           outputs= [('PATH',Path)],
           input_descriptions={'X': 'Matrix representing the data of the problem'},
           parameter_descriptions=regress_parameter_descriptions,
           output_descriptions= {
               'PATH':"Type that contains the componants of the solution related to the computation of the path, which are : 'betas,'sigmas','lambdas','method','formulation' and 'time' "
               },
           name='regress',
           description=("The function computes the constrainted_sparse_regression vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given")
           #citations=[citations['Weiss2017']]
           )


