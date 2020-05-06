import skbio 
import numpy as np

def make_lists_tree(taxa_table):

    sktree= skbio.TreeNode.from_taxonomy(taxa_table)
    #rint(sktree.ascii_art())

    POS, E, LAB = [], [], []

    def add(x, y, node, dx : float ):
        current_index = len(POS)
        POS.append([x,y])
        LAB.append(node.name)
        CHILDREN = node.children
        n = len(CHILDREN)
        if (n>0) : 
            step = dx / n
            x_new = x - dx - step
            for i in range(n)  : 
                x_new = x_new + 2*step
                E.append([current_index, len(POS)])
                add(x_new,y+1,CHILDREN[i], step)
        
    add(0.,0.,sktree, 10.)
    for i in range(len(LAB)):
        if LAB[i] is None : LAB[i] = 'None'
        if LAB[i][-1] == '_' : remove(POS,E,LAB,i)


    return np.array(POS), E, np.array(LAB)



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


