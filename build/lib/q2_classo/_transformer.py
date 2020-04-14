import numpy as np 
import zarr

def classo_to_dir(problem):
    result = CLASSOProblemDirectoryFormat()
    
    # only do it for PATH for now
    solPATH_classo, solPATH_file = problem.solution.PATH , result.PATH
    if type(solPATH_classo) != str : 
        if type(solPATH_classo.SIGMAS)==str : 
            data = zarr.array(np.array([solPATH_classo.LAMBDAS,solPATH_classo.BETAS]) )
        else :                                
            data = zarr.array(np.array([solPATH_classo.LAMBDAS,solPATH_classo.BETAS,solPATH_classo.SIGMAS]) ) 
        
        with solPATH_file.open() as fh : 
            data.write(fh,format=SolutionFormat) # is it how it is supposed to be done ??
  
    return result
