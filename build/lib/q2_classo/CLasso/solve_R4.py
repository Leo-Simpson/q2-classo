import numpy as np
import numpy.linalg as LA
from .solve_R3 import problem_R3, Classo_R3
    
'''    
Problem    :   min h_rho((Ab - y)/sigma)sigma + simga + lambda ||b||1 with C.b= 0, sigma>0 

Dimensions :   A : m*d  ;  y : m  ;  b : d   ; C : k*d

The first function compute a solution of a Lasso problem for a given lambda. The parameters are lam (lambda/lambdamax, \in [0,1]) and pb, which has to be a 'problem_LS type', which is defined bellow in order to contain all the important parameters of the problem. One can initialise it this way : pb = class_of_problem.problem(data=(A,C,y),type_of_algo). We solve the problem without normalizing anything. 
'''    


def Classo_R4(pb,lam,e=1.):
    pb_type = pb.type      # can be 'Path-Alg' or 'DR'
    (m,d,k),(A,C,y)  = pb.dim,pb.matrix
    lamb,rho  = lam * pb.lambdamax, pb.rho
    regpath = pb.regpath
    # Only alternative to 2prox : one can use the other formulation of the problem which shows that we can augment the data and then simply solve a concomitant problem
    # (we do that with the method ODE for example becasue it is pretty efficient). Unfortunately we can do that only for fixed lambda and not for any path algorithms
    # because the augmentation of the data required depends on lambda.
    if pb_type=='Path-Alg' and not regpath:
        matrix_aug = (np.concatenate((A,lamb/(2*rho)*np.eye(m)),axis=1),np.concatenate((C,np.zeros((k,m))),axis=1),y)
        pb_aug = problem_R3(matrix_aug, 'Path-Alg', e=e)
        beta,s = Classo_R3(pb_aug, lamb / pb_aug.lambdamax)
        s = s / np.sqrt(e)
        beta= beta[:d]
        return beta,s

    # Else, we do simply doulgas rachford. Hence, the prox is not so easy to compute because there is a root of polynomial of degree 3 to compute.
    # We do that in the function prox_phi_2 which use the function prox_phi_i (prox of one componant), and it uses calc_Newton which uses newton's method with good initialization.

    if (not regpath): pb.compute_param()
    proj_sigm, QA, Q1, Q2, Proj, Anorm = pb.proj_sigm, pb.QA, pb.Q1, pb.Q2, pb.Proj, pb.Anorm
    tol     = pb.tol * LA.norm(y) # tolerance rescaled
    gamma = LA.norm(y) * pb.gam / (Anorm ** 2)
    w,zerod = lamb *gamma*pb.weights, np.zeros(d) # two vectors usefull to compute the prox of f(b)= sum(wi |bi|)
    mu, c   = pb.mu, pb.c
    root    = [0.]*len(y)
    xs,nu,o,xbar,x = pb.init


    #2prox
    if (pb_type == 'DR'):
        for i in range(pb.N):
            nv_b, nv_s = x + Q1.dot(o) - QA.dot(x) - Q2.dot(x-xbar), (xs+nu)/2
            if (i>0 and LA.norm(b-nv_b)*Anorm +LA.norm(s-nv_s)<2*tol):
                if (regpath):
                            return(b,(xs,nu,o,xbar,x),sum(s)/len(s)/pb.sigmax)
                else :      return(b,sum(s)/len(s))
                        
            s,b = nv_s, nv_b
            Ab = A.dot(b)
            p1,p2,root = prox_phi_2(xs,2*Ab-o-y,gamma/c,root,rho)
            sup = [ proj_sigm(nu)-s , p1-s , p2 + y - Ab , prox(2*b-xbar,w,zerod)-b, Proj.dot(2*b-x)-b]
            xs,nu,o,xbar,x = xs+mu*sup[0] ,  nu+mu*sup[1] ,  o+mu*sup[2] ,  xbar+mu*sup[3] ,  x+mu*sup[4]
            if (LA.norm(b)+LA.norm(s)>1e6): 
                print('DIVERGENCE')
                return(b,np.sqrt(m)*sum(s)/len(s))
        print('NO CONVERGENCE')
        return(b,sum(s)/len(s))

    print('none of the cases ! ')        
    
    

    
'''
This function compute the the solution for a given path of lam : by calling the function 'algo' for each lambda with warm start, or wuth the method ODE, by computing the whole path thanks to the ODE that rules Beta and the subgradient s, and then to evaluate it in the given finite path.  
'''
    
def pathlasso_R4(pb,path,n_active=False):
    n = pb.dim[0]
    BETA,SIGMA,tol = [],[],pb.tol

    save_init = pb.init   
    pb.regpath = True
    pb.compute_param()
    for lam in path:
        X = Classo_R4(pb,lam)
        BETA.append(X[0]), SIGMA.append(X[2])
        pb.init = X[1]
        if (type(n_active)==int) : n_act = n_active
        else : n_act = n
        if(sum([ (abs(X[0][i])>1e-1) for i in range(len(X[0])) ])>=n_act):
                pb.init, BETA, SIGMA = save_init, BETA + [BETA[-1]]*(len(path)-len(BETA)),SIGMA + [SIGMA[-1]]*(len(path)-len(SIGMA))
                return(BETA,SIGMA)
            
    pb.init = save_init
    pb.regpath = False
    return(BETA,SIGMA)








'''
Class of problem : we define a type, which will contain as keys, all the parameters we need for a given problem.
'''


class problem_R4 :
    
    def __init__(self,data,algo,rho,e=1.):
        self.N = 500000
        
        (A,C,y), self.dim = data, (data[0].shape[0],data[0].shape[1],data[1].shape[0])
        self.matrix = (A,C,y)
        
        (m,d,k) = self.dim
        self.weights = np.ones(d)
        self.tol = 1e-3

        self.regpath = False
        self.name = algo + ' Concomitant Huber'
        self.type = algo          # type of algorithm used
        rho_max = LA.norm(y,np.inf)
        self.rho = rho
        self.mu  = 1.95
         
        self.c = (d/LA.norm(A,2))**2  # parameter for Concomitant problem : the matrix is scaled as c*A^2 
        self.gam = np.sqrt(d)
        e = m
        sigmax = find_sigmax(y,rho,e)
        #sigmax = LA.norm(y)/np.sqrt(m)
        self.sigmax = sigmax
        self.lambdamax = 2 * LA.norm(A.T.dot(y), np.infty) / LA.norm(y) * np.sqrt(e)
        self.lambdamax = 2*LA.norm((A.T).dot(h_prime(y/sigmax,rho)),np.infty)
        self.init = sigmax*np.ones(m),sigmax*np.ones(m),np.zeros(m), np.zeros(d), np.zeros(d)

    def compute_param(self):
        (A, C, y) = self.matrix
        m, d, k = self.dim
        self.Anorm = LA.norm(A, 'fro')
        self.Proj = proj_c(C, d)  # Proj = I - C^t . (C . C^t )^-1 . C
        self.Q1, self.Q2 = QQ(self.c, A)
        self.QA = self.Q1.dot(A)
        self.proj_sigm = lambda vect: ([max(0,sum(vect))/len(vect)]*len(vect)) # here,
                                # compared to the Muller&Combettes paper, there is a projection more on sigma =0,
                                # for numerical issues...
        



'''
Functions used in the algorithms, modules needed : 
import numpy as np
import numpy.linalg as LA
from .../class_of_problem import problem
'''


# compute the prox of the function : f(b)= sum (wi * |bi| )
def prox(b,w,zeros): return(np.minimum(b+w,zeros)+np.maximum(b-w,zeros)) 

# Compute I - C^t (C.C^t)^-1 . C : the projection on Ker(C)
def proj_c(M,d):
    if (LA.matrix_rank(M)==0):  return(np.eye(d))
    return(np.eye(d)-LA.multi_dot([M.T,np.linalg.inv(M.dot(M.T) ),M]) )



def QQ(coef,A): return(coef*(A.T).dot(LA.inv(2*np.eye(A.shape[0])+coef*A.dot(A.T))),LA.inv(2*np.eye(A.shape[1])+coef*(A.T).dot(A)))    



# Compute the real positive root of a polynomial of degree 3 in the form : X^3 + a*X - b with Newton method and a warm start (for Comcomitant problem)
def calc_Newton(a,b,root):
    er = root**3 + a*root-b
    i=0
    bound = min(np.cbrt(b), b/a)
    while (abs(er)>1e-6 and i < 20):
        if(root<0.): root, er =0., -b
        elif(root>bound):root, er =bound, bound**3 + a*bound-b
        root= root-er/(3*root**2+a)
        er = root**3 + a*root-b
        i+=1
    if(i==20):
        r = np.roots([1.,0.,a,-b])
        root = np.amax((r[np.isreal(r)]).real)
    return(root)

    
def prox_phi_2(sig,u,gamma,warm_start,rho):
    p,q, ws = np.zeros(len(u)), np.zeros(len(u)),  np.zeros(len(u))
    for i in range(len(u)):
        p[i],q[i],ws[i]=prox_phi_i(sig[i],u[i],2*gamma,warm_start[i],rho)
    return(p,q,ws)





# Each componant of the prox. Explicit formula were given in the Combettes&Muller paper.
# It is the prox of the function (u,s) --> (0.5*h_rho(u/sigma) + 0.5)*sigma
def prox_phi_i(s,u,gamma,root,rho):
    if (u==0.): return(0,0,root)
    frac = gamma*rho/abs(u)
    term = s+gamma*(rho**2-1)/2
    bool1 = (frac>=1) 
    bool2 = (abs(u)**2<=gamma*(gamma-2*s))
    bool3 = (term<=0)
    bool4 = (abs(u)>=rho*s + gamma*rho*(1+rho**2)/2)
    
    
    if   (bool1  and bool2 ): return(0.,0.,root)
    elif (bool3 and not bool1): return(0.,u*(1-frac),root) 
    elif (not bool3 and bool4):  return(term,u*(1-frac),root)
    root = calc_Newton(2*s/gamma+1,2*abs(u)/gamma,root)
    return(s+gamma*(root**2-1)/2,u-gamma*root*np.sign(u),root)


# Compute the derivative of the huber function, particulary useful for the computing of lambdamax 
def h_prime(y,rho):
    m = len(y)
    lrho = rho*np.ones(m)
    return(np.maximum(lrho,-y)+ np.minimum(y-lrho,0))



#useful for computing lambdamax.
def find_sigmax(y,rho,e):
    m,evol = len(y), True
    F = [True]*m
    if (rho > 1):
        while(evol):
            evol = False
            s = LA.norm(y[F])/np.sqrt(e-(m-sum(F))*rho**2)
            cste= rho*s
            for j in range(m):
                if (F[j] and y[j]>cste): F[j],evol= False , True
                elif (not F[j] and not y[j]>cste) : F[j],evol= True , True
        return(s)
    else:
        print('rho too little ==> sigma is always 0')




