import skbio 
import numpy as np

def make_lists_tree(taxa_table):

    sktree= skbio.TreeNode.from_taxonomy(taxa_table)
    #rint(sktree.ascii_art())

    POS, E, LAB = [], [], []

    def add(tet, r, node, dtet : float, circular = False ):
        current_index = len(POS)
        if circular : POS.append(from_polaire(r,tet))
        else : POS.append([tet,r])
        LAB.append(node.name)
        CHILDREN = node.children
        n = len(CHILDREN)
        if (n>0) : 
            step = dtet / n
            tet_new = tet - dtet - step
            for i in range(n)  : 
                tet_new = tet_new + 2*step
                E.append([current_index, len(POS)])
                add(tet_new,r+1,CHILDREN[i], step)
        
    add(0.,0.,sktree, np.pi)
    for i in range(len(LAB)):
        if LAB[i] is None : LAB[i] = 'None'
        if LAB[i][-1] == '_' : remove(POS,E,LAB,i)


    return -np.array(POS), E, np.array(LAB)



def remove(POS,E,LAB, i ):
    LAB[i] = 'None'
    POS[i] = [0.,0.]
    pred = None
    for k in range(len(E)):
        if E[k][1]==i : 
            pred = E[k][0]
            E.pop(k)
            break
    for c in E : 
        if c[0]==i : c[0] = pred


def from_polaire(r,tet):
    x = r*np.cos(tet)
    y = r*np.sin(tet)
    return [x,y]