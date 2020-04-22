import numpy as np
from q2_classo.CLasso import *
from qiime2.plugin import (SemanticType,Plugin, Int, Float, Range, Metadata, Str, Bool,
     Choices, MetadataColumn, Categorical, List,
     Citations, TypeMatch, Numeric)
from q2_types.feature_table import FeatureTable, Composition
from q2_types.feature_data import FeatureData
import qiime2



def regress(
            features : np.ndarray,
            y : qiime2.NumericMetadataColumn,
            c : np.ndarray  = None,
            #PATH parameters :
            path : bool = False,
            path_numerical_method : str         = 'not specified',
            path_n_active         : int         = 0,
            path_lambdas          : list  = None,
            path_nlam_log         : int         = 40,
            path_lamin_log        : float       = 1e-2,

            #CV parameters :
            cv : bool                         = False,
            cv_numerical_method : str         = 'not specified',
            cv_seed             : int         = None, # do something here ! for now it can be a bool !
            cv_lambdas          : list  = None, # to do 
            cv_one_se            : bool        = True,
            cv_subsets          : int         = 5,

            #StabSel parameters :
            stabsel : bool = False,
            stabsel_numerical_method : str    = 'not specified',
            stabsel_seed             : int    = None, # do something here ! for now it can be a bool !
            stabsel_lam              : float  = -1.0, # if negative, then it means 'theoretical'
            stabsel_true_lam         : bool   = True,
            stabsel_method           : str    = 'first',
            stabsel_b                : int    = 50, 
            stabsel_q                : int    = 10,
            stabsel_percent_ns       : float  = 0.5,
            stabsel_lamin            : float  = 1e-2,
            stabsel_threshold        : float  = 0.7,
            stabsel_threshold_label  : float  = 0.4, # might unneeded here, but needed for visualisation

            #LAMfixed parameters :
            lamfixed : bool = True,
            lamfixed_numerical_method : str  = 'not specified',
            lamfixed_lam              : float = -1.0, # if negative, then it means 'theoretical'
            lamfixed_true_lam         : bool  = True,
            
            #Formulation parameters
            concomitant: bool      = True,
            huber      : bool      = False,
            rho        : float     = 1.345,
            rescale    : bool      = False) -> classo_problem :


    y = y.to_series().to_numpy()

    problem = classo_problem(features, y , C = c, rescale=rescale)
    problem.formulation.huber       = huber
    problem.formulation.concomitant = concomitant
    problem.formulation.rho         = rho

    problem.model_selection.PATH= path
    if path:
        param = problem.model_selection.PATHparameters
        param.numerical_method = path_numerical_method
        param.n_active         = path_n_active
        if path_lambdas is None: 
            param.lambdas = np.array([10**(np.log10(path_lamin_log) * float(i) / path_nlam_log) for i in range(0,path_nlam_log) ] )
        else : param.lambdas=  path_lambdas

    problem.model_selection.CV = cv
    if cv :
        param = problem.model_selection.CVparameters
        param.numerical_method = cv_numerical_method
        param.seed = cv_seed    
        param.oneSE = cv_one_se  
        param.Nsubsets = cv_subsets  
        if cv_lambdas is None: param.lambdas =  np.linspace(1., 1e-3, 500)
        else                 : param.lambdas =  cv_lambdas

    problem.model_selection.StabSel = stabsel
    if stabsel : 
        param = problem.model_selection.StabSelparameters
        param.numerical_method = stabsel_numerical_method
        param.seed = stabsel_seed
        param.true_lam = stabsel_true_lam
        param.method = stabsel_method 
        param.B = stabsel_b
        param.q = stabsel_q
        param.percent_nS = stabsel_percent_ns
        param.lamin = stabsel_lamin
        param.threshold = stabsel_threshold
        param.threshold_label = stabsel_threshold_label
        if (stabsel_lam>0.): param.lam = stabsel_lam
        else               : param.lam = 'theoretical'

    problem.model_selection.LAMfixed = lamfixed
    if lamfixed: 
        param = problem.model_selection.LAMfixedparameters
        param.numerical_method = lamfixed_numerical_method
        param.true_lam = lamfixed_true_lam
        if (lamfixed_lam>0.): param.lam = lamfixed_lam
        else                : param.lam = 'theoretical'

    problem.solve()

    return problem
    
    '''
    solution_PATH, solution_CV, solution_StabSel, solution_LAM = problem.solution.PATH, problem.solution.CV, problem.solution.StabSel, problem.solution.LAMfixed

    output = dict()
    if PATH : output["PATH"] = [solution_PATH.BETAS, solution_PATH.SIGMAS, solution_PATH.LAMBDAS, solution_PATH.method, solution_PATH.formulation, solution_PATH.time]
    else : output["PATH"] = 'not_computed'

    if CV : output["CV"] = [solution_CV.beta]
    else : output["CV"] = 'not_computed'

    if StabSel : output["StabSel"] = [solution_StabSel.beta]
    else : output["StabSel"] = 'not_computed'

    if LAMfixed : output["LAMfixed"] = [solution_LAM.beta]
    else : output["LAMfixed"] = 'not_computed'
    return output
    '''

def generate_data(n : int = 100,
                  d : int = 100,
                  d_nonzero : int = 5
                    ) -> np.ndarray :

    m,d,d_nonzero,k,sigma =100,100,5,1,0.5
    (X,C,y),sol = random_data(n,d,d_nonzero,0,0.5,zerosum=True,seed= 4)
    return X






            
'''
Example from feature table :
 



def summarize(output_dir: str, table: biom.Table,
              sample_metadata: qiime2.Metadata = None) -> None:
    number_of_features, number_of_samples = table.shape

    sample_summary, sample_frequencies = _frequency_summary(
        table, axis='sample')
    if number_of_samples > 1:

        # Calculate the bin count, with a minimum of 5 bins
        IQR = sample_summary['3rd quartile'] - sample_summary['1st quartile']
        if IQR == 0.0:
            bins = 5
        else:
            # Freedman–Diaconis rule
            bin_width = (2 * IQR) / (number_of_samples ** (1/3))

            bins = max((sample_summary['Maximum frequency'] -
                        sample_summary['Minimum frequency']) / bin_width, 5)

        sample_frequencies_ax = sns.distplot(sample_frequencies, kde=False,
                                             rug=True, bins=int(round(bins)))
        sample_frequencies_ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        sample_frequencies_ax.set_xlabel('Frequency per sample')
        sample_frequencies_ax.set_ylabel('Number of samples')
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.pdf'))
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.png'))
        plt.gcf().clear()

    feature_summary, feature_frequencies = _frequency_summary(
        table, axis='observation')
    if number_of_features > 1:
        feature_frequencies_ax = sns.distplot(feature_frequencies, kde=False,
                                              rug=False)
        feature_frequencies_ax.set_xlabel('Frequency per feature')
        feature_frequencies_ax.set_ylabel('Number of features')
        feature_frequencies_ax.set_xscale('log')
        feature_frequencies_ax.set_yscale('log')
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.pdf'))
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.png'))

    sample_summary_table = q2templates.df_to_html(
        sample_summary.apply('{:,}'.format).to_frame('Frequency'))
    feature_summary_table = q2templates.df_to_html(
        feature_summary.apply('{:,}'.format).to_frame('Frequency'))

    index = os.path.join(TEMPLATES, 'summarize_assets', 'index.html')
    context = {
        'number_of_samples': number_of_samples,
        'number_of_features': number_of_features,
        'total_frequencies': int(np.sum(sample_frequencies)),
        'sample_summary_table': sample_summary_table,
        'feature_summary_table': feature_summary_table,
    }

    feature_qualitative_data = _compute_qualitative_summary(table)
    sample_frequencies.sort_values(inplace=True, ascending=False)
    feature_frequencies.sort_values(inplace=True, ascending=False)
    sample_frequencies.to_csv(
        os.path.join(output_dir, 'sample-frequency-detail.csv'))
    feature_frequencies.to_csv(
        os.path.join(output_dir, 'feature-frequency-detail.csv'))

    feature_frequencies = feature_frequencies.astype(int) \
        .apply('{:,}'.format).to_frame('Frequency')
    feature_frequencies['# of Samples Observed In'] = \
        pd.Series(feature_qualitative_data).astype(int).apply('{:,}'.format)
    feature_frequencies_table = q2templates.df_to_html(feature_frequencies)
    sample_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'sample-frequency-detail.html')
    feature_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'feature-frequency-detail.html')

    context.update({'max_count': sample_frequencies.max(),
                    'feature_frequencies_table': feature_frequencies_table,
                    'feature_qualitative_data': feature_qualitative_data,
                    'tabs': [{'url': 'index.html',
                              'title': 'Overview'},
                             {'url': 'sample-frequency-detail.html',
                              'title': 'Interactive Sample Detail'},
                             {'url': 'feature-frequency-detail.html',
                              'title': 'Feature Detail'}]})

    # Create a JSON object containing the Sample Frequencies to build the
    # table in sample-frequency-detail.html
    sample_frequencies_json = sample_frequencies.to_json()

    templates = [index, sample_frequency_template, feature_frequency_template]
    context.update({'frequencies_list':
                    json.dumps(sorted(sample_frequencies.values.tolist()))})
    if sample_metadata is not None:
        context.update({'vega_spec':
                        json.dumps(vega_spec(sample_metadata,
                                             sample_frequencies
                                             ))
                        })
    context.update({'sample_frequencies_json': sample_frequencies_json})
    q2templates.util.copy_assets(os.path.join(TEMPLATES,
                                              'summarize_assets',
                                              'vega'),
                                 output_dir)
    q2templates.render(templates, output_dir, context=context)
'''