import qiime2.plugin.model as model


class SolutionFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass

    def add(self,L):
        # method add will simply add some caracteristics to serialize the data into a text file. 
        for (name,data) in L:
            self.add_single(name,data)

    def add_single(self,name,data):
        pass

    
class CLASSOProblemDirectoryFormat(model.DirectoryFormat):
    PATH        = model.File('classo.path'       , format=SolutionFormat)
    CV          = model.File('classo.CV'         , format=SolutionFormat)
    StabSel     = model.File('classo.StabSel'    , format=SolutionFormat)
    LAMfixed    = model.File('classo.LAMfixed'   , format=SolutionFormat)

    
def classo_to_dir(problem):
    result = CLASSOProblemDirectoryFormat()



    solPATH_classo, solPATH_dir= problem.solution.PATH , result.PATH
    if type(solPATH_classo) != str : solPATH_dir.add([('betas',solPATH_classo.BETAS),
                                                      ('sigmas',solPATH_classo.SIGMAS),
                                                      ('lambdas',solPATH_classo.LAMBDAS)])

    solCV_classo, solCV_dir= problem.solution.CV , result.CV
    if type(solCV_classo) != str : solCV_dir.add([('beta',solCV_classo.beta),
                                                  ('sigma',solPATH_classo.sigma)])

    solStabSel_classo, solStabSel_dir= problem.solution.StabSel , result.StabSel
    if type(solStabSel_classo) != str : solStabSel_dir.add([('distribution',solStabSel_classo.distribution),
                                                            ('selected_parameters',solStabSel_classo.selected_param)])

    solLAMfixed_classo, solLAMfixed_dir= problem.solution.LAMfixed , result.LAMfixed
    if type(solLAMfixed_classo) != str : solLAMfixed_dir.add([('beta',solLAMfixed_classo.beta),
                                                  ('sigma',solLAMfixed_classo.sigma)])


    
    return result