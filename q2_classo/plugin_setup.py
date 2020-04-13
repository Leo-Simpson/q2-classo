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
           outputs= [('classo_problem',CLASSOProblemDirectoryFormat)],
           input_descriptions={'X': 'Matrix representing the data of the problem'},
           parameter_descriptions=regress_parameter_descriptions,
           output_descriptions= {
               'classo_problem':"Directory format that will contain all information about the problem solved"
               },
           name='regress',
           description=("The function computes the constrainted_sparse_regression vector with respect to the formulationof regression that is asked and with respect to the model selection parameters given")
           #citations=[citations['Weiss2017']]
           )



@plugin.register_transformer
def _0(obj: classo_problem) -> CLASSOProblemDirectoryFormat:
    return classo_to_dir(obj)


@plugin.register_transformer
def _1(obj: CLASSOProblemDirectoryFormat) -> classo_problem:
    return dir_to_classo(obj)