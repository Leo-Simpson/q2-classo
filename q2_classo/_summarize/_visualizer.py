import os
import json
import zarr
import pandas as pd
import numpy as np
import q2templates
import shutil
import matplotlib.pyplot as plt
from ..CLasso.misc_functions import influence, colo

dir_path = os.path.dirname(os.path.realpath(__file__))
assets = os.path.join(dir_path, 'assets')
dir_form = os.path.join(dir_path, 'form')




def summarize(output_dir: str, problem : zarr.hierarchy.Group):

    context = build_context(output_dir,problem)

    index = os.path.join(assets, 'index.html')
    overview_template = os.path.join(assets, 'overview.html')
    path_template = os.path.join(assets, 'path.html')
    cv_template = os.path.join(assets, 'cv.html')
    stabsel_template = os.path.join(assets, 'stabsel.html')
    lam_fixed_template = os.path.join(assets, 'lam-fixed.html')
    templates = [index, overview_template, path_template, cv_template, stabsel_template, lam_fixed_template]
    q2templates.render(templates, output_dir, context=context)
    


def build_context(output_dir,problem):

    labels = problem['label']
    features = pd.DataFrame(problem['data/X'],columns=labels)
    y = pd.DataFrame({'y':problem['data/y']})
    c = pd.DataFrame(problem['data/C'], columns = labels)
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
        'n' : len(problem['data/X']),
        'd' : len(problem['data/X'][0]),
        'k' : len(problem['data/C'])
    }
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

        plot_path(np.array(problem['solution/PATH/BETAS']),SIGMAS,problem['solution/PATH/LAMBDAS'],output_dir,labels)

    if dico_ms['CV']:
        context['cv'] = True
        xGraph, yGraph, standard_error = problem['solution/CV/xGraph'], problem['solution/CV/yGraph'], problem['solution/CV/standard_error']
        context['tabs'].append({'title': 'Cross-Validation','url': 'cv.html'})
        dico_cv = { **problem['model_selection/CVparameters'].attrs.asdict(), **problem['solution/CV'].attrs.asdict() }
        dico_cv['lamin']= min(xGraph)
        dico_cv['Nlam']= len(xGraph)
        beta = pd.DataFrame(data={'label':problem['label'],'beta':problem['solution/CV/refit']})
        beta.to_csv(os.path.join(output_dir,'CV-beta.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/CV/selected_param'])
        beta_support = beta[selected_param]
        dico_cv['htmlbeta']=q2templates.df_to_html(beta_support, index=False)

        context['dicocv']= dico_cv

        if (dico_cv['oneSE']): lam = dico_cv['lambda_1SE']
        else : lam = dico_cv['lambda_min']

        plot_beta(np.array(problem['solution/CV/refit']),selected_param,output_dir,labels,'cv-refit.png',r"Refitted coefficients of $\beta$ after CV model selection finds $\lambda$ = "+str(lam))
        plot_cv(xGraph, yGraph,dico_cv['index_1SE'],dico_cv['index_min'],standard_error, output_dir, 'cv-graph.png')

    if dico_ms['StabSel']:
        context['stabsel'] = True
        context['tabs'].append({'title': 'Stability Selection','url': 'stabsel.html'})
        dico_stabsel = { **problem['model_selection/StabSelparameters'].attrs.asdict(), **problem['solution/StabSel'].attrs.asdict() }
        
        stability = pd.DataFrame(data={'label':problem['label'],'stability-probability':problem['solution/StabSel/distribution']})
        stability.to_csv(os.path.join(output_dir,'StabSel-prob.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/StabSel/selected_param'])
        stability_support = stability[selected_param]

        dico_stabsel['nsel']=len(stability_support)
        dico_stabsel['htmlstab']=q2templates.df_to_html(stability_support, index=False)

        context['dicostabsel']= dico_stabsel

        plot_beta(np.array(problem['solution/StabSel/refit']),selected_param,output_dir,labels,'stabsel-refit.png',r"Refitted coefficients of $\beta$ after stability selection")
        plot_stability(problem['solution/StabSel/distribution'], selected_param, dico_stabsel['threshold'],dico_stabsel['method'], labels, output_dir,'stabsel-graph.png')

    if dico_ms['LAMfixed']:
        context['lam'] = True
        context['tabs'].append({'title': 'LAM fixed','url': 'lam-fixed.html'})
        dico_lam = { **problem['model_selection/LAMfixedparameters'].attrs.asdict(), **problem['solution/LAMfixed'].attrs.asdict() }
        dico_lam['lamtype']= problem['model_selection/LAMfixedparameters'].attrs['lam']

        beta = pd.DataFrame(data={'label':problem['label'],'beta':problem['solution/LAMfixed/refit']})
        beta.to_csv(os.path.join(output_dir,'LAM-beta.csv'),header=True, index=False)
        selected_param = np.array(problem['solution/LAMfixed/selected_param'])
        beta_support = beta[selected_param]
        dico_lam['htmlbeta']=q2templates.df_to_html(beta_support, index=False)

        context['dicolam']= dico_lam

        plot_beta(np.array(problem['solution/LAMfixed/beta']),None,output_dir,labels,'lam-beta.png',r"Coefficients of $\beta$ at $\lambda$ = "+str(dico_lam['lam']))
        plot_beta(np.array(problem['solution/LAMfixed/refit']),selected_param,output_dir,labels,'lam-refit.png',r"Reffited coefficients of $\beta$ at $\lambda$ = "+str(dico_lam['lam']))

    return context











def name_formulation(dictio,output_dir):
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


def plot_cv(xGraph, yGraph,index_1SE, index_min,SE, output_dir, name):
    fig = plt.figure()
    ax = fig.subplots()

    mse_max = 1.
    for j in range(len(xGraph)):
        if (yGraph[j] < mse_max): break

    ax.errorbar(xGraph[j:], yGraph[j:], SE[j:], label='mean over the k groups of data', errorevery = 10 )
    ax.axvline(x=xGraph[index_min], color='k', label=r'$\lambda$ (min MSE)')
    ax.axvline(x=xGraph[index_1SE],color='r',label=r'$\lambda$ (1SE) ')
    ax.set_title(r" " )
    ax.set_xlabel(r"$\lambda / \lambda_{max}$")
    ax.set_ylabel(r"Mean-Squared Error (MSE) ")
    ax.legend()

    fig.savefig(os.path.join(output_dir, name))


def plot_stability(distribution, selected_param, threshold, method, labels, output_dir,name):
    fig = plt.figure()
    ax = fig.subplots()

    D, selected = np.array(distribution), selected_param
    unselected = [not i for i in selected]
    Dselected, Dunselected  = np.zeros(len(D)), np.zeros(len(D))
    Dselected[selected], Dunselected[unselected] = D[selected], D[unselected]

    ax.bar(range(len(Dselected)), Dselected, color='r', label='selected coefficients')
    ax.bar(range(len(Dunselected)), Dunselected, color='b', label='unselected coefficients')
    ax.axhline(y=threshold, color='g',label='Threshold : thresh = '+ str(threshold))

    ax.set_xticks(np.where(selected_param)[0])
    ax.set_xticklabels(np.array(labels)[selected_param], rotation=30)
    
    ax.set_xlabel(r"Coefficient index $i$")
    ax.set_ylabel(r"Selection probability ")
    ax.set_title(r"Stability selection profile with method "+ method)
    ax.legend()


    fig.savefig(os.path.join(output_dir, name))

