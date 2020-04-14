import qiime2.plugin.model as model
from q2_types.feature_table._transformer import _4 as biom_to_pd

class SolutionFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass

class CLASSOProblemDirectoryFormat(model.DirectoryFormat):
    PATH        = model.File('classo.path'       , format=SolutionFormat)
    CV          = model.File('classo.CV'         , format=SolutionFormat)
    StabSel     = model.File('classo.StabSel'    , format=SolutionFormat)
    LAMfixed    = model.File('classo.LAMfixed'   , format=SolutionFormat)



'''
function for transformer
'''


import numpy as np 
import zarr


def write_with_numpy(file,data):
    # function to write !! 
    path = file.path_maker()
    array = zarr.array(data, dtype=float)
    zarr.save(str(path),array)


def classo_to_dir(problem):
    result = CLASSOProblemDirectoryFormat()
    
    # only do it for PATH for now
    solPATH_classo= problem.solution.PATH
    if type(solPATH_classo) != str : 
        if type(solPATH_classo.SIGMAS)==str : data = np.array([solPATH_classo.LAMBDAS,solPATH_classo.BETAS])
        else :   data = np.array([solPATH_classo.LAMBDAS,solPATH_classo.BETAS,solPATH_classo.SIGMAS])

        write_with_numpy(result.PATH,data)

    solCV_classo= problem.solution.CV
    if type(solCV_classo) != str : 
        data = solCV_classo.beta # to do here
        write_with_numpy(result.CV,data)

    solStabSel_classo= problem.solution.StabSel
    if type(solStabSel_classo) != str : 
        data = solStabSel_classo.distribution # to do here
        write_with_numpy(result.StabSel,data)

    solLAMfixed_classo= problem.solution.LAMfixed
    if type(solLAMfixed_classo) != str : 
        data = solLAMfixed_classo.beta # to do here
        write_with_numpy(result.LAMfixed,data)
        
    return result



def table_to_array(biom_obj):
    data_frame = biom_to_pd(biom_obj)
    array = np.array(data_frame)
    return array
