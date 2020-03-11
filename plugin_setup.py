from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch)


citations = Citations.load('citations.bib', package='q2_feature_table')
plugin = Plugin(
name='classo',
version=q2_feature_table.__version__,
website='https://github.com/Leo-Simpson/q2-classo',
package='q2_classo',
short_description=('Package for constrained sparse regression and classification'),
description=('This is QIIME 2 plugin that enables sparse and robust linear regression and classification with linear equality constraints on the model parameters.')
)

plugin.methods.register_function(
           function=somewhere.regress,
           inputs={'name': type, },
           parameters={'name1': type, 'name2': type },
           outputs=[('name',type)],
           input_descriptions={'name': 'Description of name'},
           parameter_descriptions={
           'name1': ('Descrption of parameter name1'),
           'name2': ('Descrption of parameter name2')
           },
           output_descriptions={
           'rarefied_table': 'The resulting rarefied feature table.'
           },
           name='Rarefy table',
           description=("Description of the function"),
           citations=[citations['Weiss2017']]
           )

plugin.methods.register_function(
                                 function=somewhere.classify,
                                 inputs={'name': type, },
                                 parameters={'name1': type, 'name2': type },
                                 outputs=[('name',type)],
                                 input_descriptions={'name': 'Description of name'},
                                 parameter_descriptions={
                                 'name1': ('Descrption of parameter name1'),
                                 'name2': ('Descrption of parameter name2')
                                 },
                                 output_descriptions={
                                 'rarefied_table': 'The resulting rarefied feature table.'
                                 },
                                 name='Rarefy table',
                                 description=("Description of the function"),
                                 citations=[citations['Weiss2017']]
                                 )

