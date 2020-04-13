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



    parameters = model.File('parameters'    , format= ParametersFormat)
    solution  = model.File('solution'     , format= SolutionFormat)

    
