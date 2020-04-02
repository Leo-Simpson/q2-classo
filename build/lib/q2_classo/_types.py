from qiime2.plugin import SemanticType
import qiime2.plugin.model as model

from q2_classo.plugin_setup import plugin

CLASSOProblem = SemanticType("CLASSOProblem")

ConstraintMatrix = SemanticType("ConstraintMatrix")

plugin.register_semantic_types(CLASSOProblem)



class CLASSOProblemDirectoryFormat(model.DirectoryFormat):
    pass # to do 


plugin.register_semantic_type_to_format(CLASSOProblem, 
                                        artifact_format=CLASSOProblemDirectoryFormat)
