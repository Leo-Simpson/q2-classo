import skbio 
import numpy as np

def make_lists_tree(sktree):

    #sktree= skbio.TreeNode.from_taxonomy(taxa_table)
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



def build_subtree(sktree,label_node): 
    sub_sktree = sktree  # shir ? http://scikit-bio.org/docs/0.5.1/generated/skbio.tree.TreeNode.shear.html#skbio.tree.TreeNode.shear
    return(sub_sktree)



def tree_to_matrix(tree,label, with_repr = False):
    # to do here : given a skbio tree and the beta-labels in a given order, return the matrix A and the new labels corresponding
    dicti = dict()
    d = len(label)
    LEAVES = [tip.name for tip in tree.tips()]
    order = [] # list that will give the order in which the codes are added in the dicti, such that it will be easy to remove similar nodes
    for i in range(d): 
        name_leaf = label[i]
        dicti[name_leaf] = np.zeros(d)
        dicti[name_leaf][i] = 1
        order.append(name_leaf)
        if name_leaf in LEAVES :
            for n in tree.find(name_leaf).ancestors() : 
                ancest = n.name
                if ancest[-1] != '_': 
                    if not ancest in dicti : 
                        dicti[ancest] = np.zeros(d)
                        order.append(ancest)
                    dicti[ancest][i] = 1


    L,label2 = [], []


    for node in tree.levelorder():
        nam = node.name
        if nam in dicti and not nam in label2 : 
            label2.append(nam)
            L.append(dicti[nam])

    to_keep = remove_same_vect(L , label2,order)


    return np.array(L)[to_keep].T, np.array(label2)[to_keep]


def remove_same_vect(L , label, order): 
    K = len(L)
    to_keep = np.array([True]*K)
    j = label.index(order[0] )
    col = L[j]
    for i in range(K-1):
        new_j = label.index(order[i+1])
        new_col = L[new_j]
        if np.array_equal(col,new_col) : 
            to_keep[new_j]=False
        else : j, col = new_j, new_col

    return to_keep