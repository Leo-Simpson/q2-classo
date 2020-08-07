import os
import json
import zarr
import pandas as pd
import numpy as np
import q2templates
import shutil
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from plotly import graph_objects, express, offline
import skbio

from time import time
from math import ceil


from .._tree import make_lists_tree, plot_tree

dir_path = os.path.dirname(os.path.realpath(__file__))
assets = os.path.join(dir_path, 'assets')
dir_form = os.path.join(dir_path, 'form')

colors = {"threshold":"red", "selected":"green", "unselected":"blue"}


def summarize(output_dir: str, problem : zarr.hierarchy.Group, taxa : skbio.TreeNode = None, maxplot : int =200):
    context = build_context(output_dir,problem, taxa,maxplot)

    index = os.path.join(assets, 'index.html')
    overview_template = os.path.join(assets, 'overview.html')
    path_template = os.path.join(assets, 'path.html')
    cv_template = os.path.join(assets, 'cv.html')
    stabsel_template = os.path.join(assets, 'stabsel.html')
    lam_fixed_template = os.path.join(assets, 'lam-fixed.html')
    
    templates = [index, overview_template, path_template, cv_template, stabsel_template, lam_fixed_template]
    
    q2templates.render(templates, output_dir, context=context)
    


def build_context(output_dir,problem,taxa,max_number):
    t = time()


    labels = np.array(problem['data/label'])
    if problem['formulation'].attrs['intercept'] : 
        features, c = pd.DataFrame(problem['data/X'],columns=labels[1:]), pd.DataFrame(problem['data/C'],columns=labels[1:]) 
    else :
        features, c = pd.DataFrame(problem['data/X'],columns=labels), pd.DataFrame(problem['data/C'],columns=labels) 
    y = pd.DataFrame({'y':problem['data/y']})

    #tree = problem['data'].attrs['tree']
    

    features.to_csv(os.path.join(output_dir, 'features.csv'),header=True, index=False)
    y.to_csv(os.path.join(output_dir, 'samples.csv'),header=True, index=False)
    c.to_csv(os.path.join(output_dir, 'constraints.csv'),header=True, index=False)


    context = {
        'path': False,'cv': False,'stabsel':False,'lam' : False,
        'labels':labels,
        'tabs': [{'title': 'Overview','url': 'overview.html'}]
    }
    dico = {
        'formulation':name_formulation(problem['formulation'].attrs.asdict(),output_dir),
        'concomitant': problem['formulation'].attrs['concomitant'],
        'with_tree' : not taxa is None,
        #'tree'  : tree,
        'n' : len(problem['data/X']),
        'd' : len(problem['data/X'][0]),
        'k' : len(problem['data/C'])
    }

    
    if dico['with_tree'] : 
        fig = plot_tree(taxa, labels)
        offline.plot(fig, filename = os.path.join(output_dir, 'tree.html'), auto_open=False)                               
        

    context['dico']=dico
    dico_ms = problem['model_selection'].attrs.asdict()


    if dico_ms['PATH']:
        context['path'] = True
        context['tabs'].append({'title': 'Lambda-path','url': 'path.html'})
        dico_path = { **problem['model_selection/PATHparameters'].attrs.asdict(), **problem['solution/PATH'].attrs.asdict() }
        dico_path['lambdas']= problem['solution/PATH/LAMBDAS']
        dico_path['lamin']= min(dico_path['lambdas'])
        dico_path['Nlam'] = len(dico_path['lambdas'])
        data = pd.DataFrame(np.array(problem['solution/PATH/BETAS']),index=dico_path['lambdas'], columns=labels)
        data.to_csv(os.path.join(output_dir, 'path.csv'),header=True, index=True)
        SIGMAS = None
        if dico['concomitant']: 
            SIGMAS = problem['solution/PATH/SIGMAS']

        context['dicopath']= dico_path

        plot_path(np.array(problem['solution/PATH/BETAS']),SIGMAS,problem['solution/PATH/LAMBDAS'],output_dir,labels,"beta-path.html","sigma-path.html")


    if dico_ms['CV']:
        context['cv'] = True
        xGraph, yGraph, standard_error = problem['solution/CV/xGraph'], problem['solution/CV/yGraph'], problem['solution/CV/standard_error']
        context['tabs'].append({'title': 'Cross-Validation','url': 'cv.html'})
        dico_cv = { **problem['model_selection/CVparameters'].attrs.asdict(), **problem['solution/CV'].attrs.asdict() }
        dico_cv['lamin']= min(xGraph)
        dico_cv['Nlam']= len(xGraph)
        beta = pd.DataFrame(data={'label':labels,'beta':problem['solution/CV/refit']})
        beta.to_csv(os.path.join(output_dir,'CV-beta.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/CV/selected_param'])
        beta_support = beta[selected_param]
        dico_cv['htmlbeta']=q2templates.df_to_html(beta_support, index=False)

        context['dicocv']= dico_cv

        if (dico_cv['oneSE']): lam = dico_cv['lambda_1SE']
        else : lam = dico_cv['lambda_min']

        plot_beta(np.array(problem['solution/CV/refit']),output_dir,labels,'cv-refit.html',r"Refitted coefficients of beta after CV model selection finds lambda = "+str(lam),max_number)
        plot_cv(xGraph, yGraph,dico_cv['index_1SE'],dico_cv['index_min'],standard_error, output_dir, 'cv-graph.html')


    if dico_ms['StabSel']:
        context['stabsel'] = True
        context['tabs'].append({'title': 'Stability Selection','url': 'stabsel.html'})
        dico_stabsel = { **problem['model_selection/StabSelparameters'].attrs.asdict(), **problem['solution/StabSel'].attrs.asdict() }
        dico_stabsel["with_path"] = ( dico_stabsel['method'] == 'first'  )

        stability = pd.DataFrame(data={'label':labels,'stability-probability':problem['solution/StabSel/distribution']})
        stability.to_csv(os.path.join(output_dir,'StabSel-prob.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/StabSel/selected_param'])
        stability_support = stability[selected_param]
        
        if dico['with_tree'] : 
            fig = plot_tree(taxa, labels, selected_labels = labels[selected_param])
            offline.plot(fig, filename = os.path.join(output_dir, 'StabSel-tree.html'), auto_open=False)  

        dico_stabsel['nsel']=len(stability_support)
        dico_stabsel['htmlstab']=q2templates.df_to_html(stability_support, index=False)

        context['dicostabsel']= dico_stabsel
        plot_beta(np.array(problem['solution/StabSel/refit']),output_dir,labels,'stabsel-refit.html',r"Refitted coefficients of beta after stability selection",max_number)
        plot_stability(np.array(stability['stability-probability']), selected_param, dico_stabsel['threshold'],dico_stabsel['method'], labels, output_dir,'stabsel-graph.html', max_number)
        
        if dico_stabsel["with_path"]: 
            plot_stability_path(problem['solution/StabSel/lambdas_path'], problem['solution/StabSel/distribution_path'],
                                selected_param,dico_stabsel['threshold'],dico_stabsel['method'],labels, output_dir, 'stabsel-path.html')



    if dico_ms['LAMfixed']:
        context['lam'] = True
        context['tabs'].append({'title': 'LAM fixed','url': 'lam-fixed.html'})
        dico_lam = { **problem['model_selection/LAMfixedparameters'].attrs.asdict(), **problem['solution/LAMfixed'].attrs.asdict() }
        dico_lam['lamtype']= problem['model_selection/LAMfixedparameters'].attrs['lam']

        beta = pd.DataFrame(data={'label':labels,'beta':problem['solution/LAMfixed/refit']})
        beta.to_csv(os.path.join(output_dir,'LAM-beta.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/LAMfixed/selected_param'])
        beta_support = beta[selected_param]
        dico_lam['htmlbeta']=q2templates.df_to_html(beta_support, index=False)

        context['dicolam']= dico_lam


        plot_beta(np.array(problem['solution/LAMfixed/beta']),output_dir,labels,'lam-beta.html',r"Coefficients of beta at lambda = "+str(dico_lam['lam']),max_number)
        plot_beta(np.array(problem['solution/LAMfixed/refit']),output_dir,labels,'lam-refit.html',r"Reffited coefficients of beta at lambda = "+str(dico_lam['lam']),max_number)



    return context











def name_formulation(dictio,output_dir):
    if dictio['classification']:
        if dictio['concomitant']:
            shutil.copy(os.path.join(dir_form, 'C2.png'),os.path.join(output_dir, 'formula.png'))
            return 'C2 (classification and huber with rho = '+str(dictio['rho_classification']) +' )'
        else : 
            shutil.copy(os.path.join(dir_form, 'C1.png'),os.path.join(output_dir, 'formula.png'))
            return 'C1 (classification)'
        
    else:
        if dictio['concomitant']:
            if dictio['huber']:
                shutil.copy(os.path.join(dir_form, 'R4.png'),os.path.join(output_dir, 'formula.png'))
                return 'R4 (concomitant and huber with e = '+ str(dictio['e']) +' and rho = '+str(dictio['rho']) +' )'
            else :
                shutil.copy(os.path.join(dir_form, 'R3.png'),os.path.join(output_dir, 'formula.png'))
                return 'R3 (concomitant with e = ' + str(dictio['e']) + ' )'
        else : 
            if dictio['huber']:
                shutil.copy(os.path.join(dir_form, 'R2.png'),os.path.join(output_dir, 'formula.png'))
                return 'R2 (huber with rho = '+str(dictio['rho']) +' )'
            else :
                shutil.copy(os.path.join(dir_form, 'R1.png'),os.path.join(output_dir, 'formula.png'))
                return 'R1 (classic lasso formulation)'





def plot_path(BETAS, SIGMAS, LAMBDAS, directory, labels, name1, name2):
    fig = graph_objects.Figure(layout_title_text=r"Coefficients across lambda-path")
    maxB = np.amax(abs(BETAS))
    start = int(labels[0]=="intercept")
    for i in range(start,len(BETAS[0])):
        if np.amax(abs(BETAS[:,i]))>maxB*0.01 : 
            fig.add_trace(graph_objects.Scatter(x=LAMBDAS, y=BETAS[:,i],
                                name=labels[i]))
    fig.update_xaxes(title_text=r"lambda")
    fig.update_yaxes(title_text=r"Coefficients beta_i ")
    offline.plot(fig, filename = os.path.join(directory, name1), auto_open=False)


    if not SIGMAS is None:
        fig2 = graph_objects.Figure(layout_title_text=r"Scale estimate across lambda-path")
        fig2.add_trace(graph_objects.Scatter(x=LAMBDAS, y=SIGMAS,
                                name="sigma"))
        fig2.update_xaxes(title_text=r"lambda")
        fig2.update_yaxes(title_text=r"Scale sigma ")
        offline.plot(fig2, filename = os.path.join(directory, name2), auto_open=False)





def plot_beta(beta,directory,labels,name,title,max_number):

    reduc = np.ones(len(beta), dtype='bool')
    if sum(beta != 0.) > max_number : 
        reduc[np.argpartition(abs(beta), -max_number)[:-max_number] ] = False
    else :
        reduc[beta==0.] = False
    if labels[0]=="intercept": reduc[0]=False

    data = {'index': range(len(beta[reduc])), "Coefficient i of beta": beta[reduc] , 'label': labels[reduc] }
    fig = express.bar(data, x='index', y="Coefficient i of beta", hover_data=['label'])
    fig.update_layout(title= title)

    offline.plot(fig, filename = os.path.join(directory, name), auto_open=False)


def plot_cv(xGraph, yGraph,index_1SE, index_min,SE, directory, name):
    mse_max, j = 10*SE[index_min], 0
    jmax = len(yGraph)-1
    while(yGraph[jmax]>100*yGraph[index_min]) : jmax-=1
    while ( j < index_1SE - 30 and yGraph[j] > mse_max) : j+=1

    y_max = max(yGraph[j:jmax])
    fig = graph_objects.Figure()
    
    #errors = np.zeros(len(SE[j:]))
    #errors[np.arange(0,len(errors),ceil(len(errors)/20))] = SE[j:][np.arange(0,len(errors),ceil(len(errors)/20))]
    errors = SE[j:]

    fig.add_trace(graph_objects.Scatter(x=xGraph[j:jmax], y=yGraph[j:jmax], name = "MSE",
                                error_y=dict(
                                type='data', # value of error bar given in data coordinates
                                array=errors,
                                visible=True)
                                ))
    
    fig.add_trace(graph_objects.Scatter(x=[xGraph[index_min],xGraph[index_min]], y=[0,y_max], mode = "lines",
                                name="Lambda min MSE"))
    fig.add_trace(graph_objects.Scatter(x=[xGraph[index_1SE],xGraph[index_1SE]], y=[0,y_max], mode = "lines",
                                name="Lambda 1SE"))
    
    fig.update_xaxes(title_text="lambda / lambda_max ")
    fig.update_yaxes(title_text="Mean-Squared Error (MSE) ")

    offline.plot(fig, filename = os.path.join(directory, name), auto_open=False)    



def plot_stability(distribution, selected_param, threshold, method, labels, directory, name, max_number):
    reduc = np.ones(len(distribution), dtype='bool')
    if len(distribution) > max_number : 
        reduc[np.argpartition(distribution, -max_number)[:-max_number] ] = False
    if labels[0]=="intercept": reduc[0]=False


    data = {'index': range(len(distribution[reduc])), "Selection probability": distribution[reduc] , 'label': labels[reduc], 'selected':selected_param[reduc] }
    fig = express.bar(data, x='index', y="Selection probability", hover_data=['selected','label'], color='selected',color_discrete_map={True:colors["selected"],False:colors["unselected"]})
    fig.update_layout(shapes=[
                    dict(type= 'line', y0= threshold, y1= threshold, x0= 0, x1= len(distribution[reduc]),line_color=colors["threshold"] )  
                            ])
    offline.plot(fig, filename = os.path.join(directory, name), auto_open=False)


def plot_stability_path(lambdas, D_path, selected, threshold,method,labels,directory,name):
    N = len(lambdas)
    D_path = np.array(D_path)
    reduc = D_path[-1,:]>0.05
    if labels[0]=="intercept": reduc[0]=False
    reduc_selected  = selected[reduc]
    d = len(reduc_selected)
    reduc_D_path = D_path[:,reduc]
    data = { "lambda":np.repeat(lambdas,d), "Selection Probability":[],"selected":list(reduc_selected)*N,"labels":list(labels[reduc])*N}
    for i in range(N):
        data["Selection Probability"].extend(reduc_D_path[i])
    fig = express.line(data, x = "lambda", y = "Selection Probability",color="selected", line_group="labels",hover_name="labels",color_discrete_map={True:colors["selected"],False:colors["unselected"]})
    fig.update_layout(title = "Stability selection profile across lambda-path with method "+ method, 
                    shapes=[ dict(type= 'line', y0= threshold, y1= threshold, x0= 0, x1= 1,line_color=colors["threshold"] ) ] 
                            )
    offline.plot(fig, filename = os.path.join(directory, name), auto_open=False)




'''
def plot_path(BETAS, SIGMAS, LAMBDAS,output_dir, labels):
    fig = plt.figure()
    ax = fig.subplots()
    l_index = influence(BETAS, 10)
    j=0
    for i in range(len(BETAS[0])):
        if j<len(l_index) and i==l_index[j]: 
            label = str(labels[i])
            j+=1
        else : 
            label=None
        ax.plot(LAMBDAS,BETAS[:,i],label=label,color=colo[i])
    ax.set_title(r"Coefficients across $\lambda$-path")
    ax.legend(loc=4, borderaxespad=0.)
    ax.set_xlabel(r"$\lambda$"), ax.set_ylabel(r"Coefficients $\beta_i$ ")
    
    fig.savefig(os.path.join(output_dir, 'beta-path.png'))
    if not SIGMAS is None:
        fig2 = plt.figure()
        ax2 = fig2.subplots()
        ax2.plot(LAMBDAS, SIGMAS, color='blue')
        ax2.set_title(r"Scale estimate across $\lambda$-path"),ax2.set_xlabel(r"$\lambda$"), ax2.set_ylabel(r"Scale $\sigma$ ")
        fig2.savefig(os.path.join(output_dir, 'sigma-path.png'))


def plot_cv(xGraph, yGraph,index_1SE, index_min,SE, output_dir, name):
    fig = plt.figure()
    ax = fig.subplots()

    mse_max = 10*SE[index_min]
    j = 0
    while ( j < index_1SE - 30 and yGraph[j] > mse_max) : j+=1

    ax.errorbar(xGraph[j:], yGraph[j:], SE[j:], label='mean over the k groups of data', errorevery = 10 )
    ax.axvline(x=xGraph[index_min], color='k', label=r'$\lambda$ (min MSE)')
    ax.axvline(x=xGraph[index_1SE],color='r',label=r'$\lambda$ (1SE) ')
    ax.set_title(r" " )
    ax.set_xlabel(r"$\lambda / \lambda_{max}$")
    ax.set_ylabel(r"Mean-Squared Error (MSE) ")
    ax.legend()

    fig.savefig(os.path.join(output_dir, name))


def plot_beta(beta,selected_param,output_dir,labels,name,title):
    fig = plt.figure()
    ax = fig.subplots()
    ax.bar(range(len(beta)), beta)
    ax.set_title(title)
    ax.set_xlabel(r"Feature of index $i$" )
    ax.set_ylabel(r"Coefficients $\beta_i$ ")
    ax.axhline(xmax=len(beta),color='k')
    if not selected_param is None : 
        ax.set_xticks(np.where(selected_param)[0])
        ax.set_xticklabels(np.array(labels)[selected_param], rotation=30)
    fig.savefig(os.path.join(output_dir, name))



def plot_stability(distribution, selected_param, threshold, method, labels, output_dir,name):
    fig = plt.figure()
    ax = fig.subplots()

    D, selected = np.array(distribution), np.array(selected_param)
    unselected = [not i for i in selected]
    Dselected, Dunselected  = np.zeros(len(D)), np.zeros(len(D))
    Dselected[selected], Dunselected[unselected] = D[selected], D[unselected]

    ax.bar(range(len(Dselected)), Dselected, color='r', label='selected coefficients')
    ax.bar(range(len(Dunselected)), Dunselected, color='b', label='unselected coefficients')
    ax.axhline(y=threshold, color='g',label='Threshold : thresh = '+ str(threshold))

    ax.set_xticks(np.where(selected_param)[0])
    ax.set_xticklabels(labels[selected_param], rotation=30)
    
    ax.set_xlabel(r"Coefficient index $i$")
    ax.set_ylabel(r"Selection probability ")
    ax.set_title(r"Stability selection profile with method "+ method)
    ax.legend()


    fig.savefig(os.path.join(output_dir, name))



def plot_stability_path(lambdas, D_path, selected, threshold,method,output_dir,name):
    fig = plt.figure()
    ax = fig.subplots()
    N = len(D_path)
    for i in range(len(selected)):
        if selected[i]: c='red'
        else          : c='blue'
        ax.plot(lambdas, [D_path[j][i] for j in range(N)], c)
    ax.axhline(y=threshold,color='green')

    p1 = mpatches.Patch(color='red', label='selected coefficients')
    p2 = mpatches.Patch(color='blue',label='unselected coefficients')
    p3 = mpatches.Patch(color='green',label='Threshold : thresh = '+ str(threshold))
    ax.legend(handles=[p1, p2, p3], loc=1)
    

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"Selection probability ")
    ax.set_title(r"Stability selection profile across $\lambda$-path with method "+ method)


    fig.savefig(os.path.join(output_dir, name))



'''


