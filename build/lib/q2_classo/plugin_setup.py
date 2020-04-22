from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric, SemanticType)
from q2_types.feature_table import FeatureTable, Composition, BIOMV210Format
from . import  *
import numpy as np
import biom
import zarr

#citations = Citations.load('citations.bib', package='q2_classo') 
plugin = Plugin(
name='classo',
                version='0.0.0.dev0',
website='https://github.com/Leo-Simpson/q2-classo',
package='q2-classo',
short_description=('Package for constrained sparse regression and classification'),
description=('This is QIIME 2 plugin that enables sparse and robust linear regression and classification with linear equality constraints on the model parameters.')
)


CLASSOProblem    = SemanticType("CLASSOProblem")
ConstraintMatrix = SemanticType("ConstraintMatrix")


plugin.register_formats(CLASSOProblemDirectoryFormat)
plugin.register_semantic_types(ConstraintMatrix)
plugin.register_semantic_type_to_format(CLASSOProblem, 
                                        artifact_format=CLASSOProblemDirectoryFormat)

@plugin.register_transformer
def _0(obj: classo_problem) -> CLASSOProblemDirectoryFormat :
    ff = CLASSOProblemDirectoryFormat()
    store = zarr.DirectoryStore(str(ff.path/'problem'))
    root = zarr.group(store=store)
    to_zarr(obj,'problem',root)
    print(root.tree())
    return ff 


@plugin.register_transformer
def _1(ff: BIOMV210Format) -> np.ndarray:
    with ff.open() as fh:
        table = biom.Table.from_hdf5(fh)
    array = table.matrix_data.toarray().T # numpy array
    return array

def _2(obj : ZarrProblemFormat) -> classo_problem : 
    store = zarr.ZipStore(obj, mode='r') #+'.zip' ???  
    root = zarr.group(store=store)
    return zarr_to_classo(root)

def _3(obj : ZarrProblemFormat) -> zarr.hierarchy.Group : 
    store = zarr.ZipStore(str(obj), mode='r') #+'.zip' ???  
    return zarr.group(store=store)

def _4(obj : ConstraintMatrix) -> np.ndarray : 
    pass



plugin.methods.register_function(
           function=regress,
           inputs={'features': FeatureTable[Composition], 'c':ConstraintMatrix},
           parameters=regress_parameters,
           outputs= [('result',CLASSOProblem)],
           input_descriptions={'features': 'Matrix representing the data of the problem','c': 'Constraint matrix, default is the zero-sum',},
           parameter_descriptions=regress_parameter_descriptions,
           output_descriptions= {
               'result':"Directory format that will contain all information about the problem solved"
               },
           name='regress',
           description=("The function computes the constrainted_sparse_regression vector with respect to the formulation of regression that is asked and with respect to the model selection parameters given")
           #citations=[citations['Weiss2017']]
           )













plugin.methods.register_function(
           function=generate_data,
           inputs={},
           parameters={'n':Int, 'd':Int, 'd_nonzero':Int},
           outputs= [('x',FeatureTable[Composition])],
           input_descriptions={},
           parameter_descriptions={'n': 'number of sample', 'd': 'number of features','d_nonzero': 'number of non nul componants in beta' },
           output_descriptions= {'x': 'Matrix representing the data of the problem'},
           name='generate_data',
           description=("Function that build random data")
           )



