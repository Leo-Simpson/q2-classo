import os
import json
import zarr
import pandas as pd
import numpy as np
import q2templates
import pkg_resources

TEMPLATES = pkg_resources.resource_filename('q2_classo','_summarize')



def summarize(output_dir: str, problem : zarr.hierarchy.Group):
    print(TEMPLATES)

    beta = pd.DataFrame(data={'label':problem['label'],'beta':problem['solution/LAMfixed/refit']})
    beta.to_csv(os.path.join(output_dir,'beta.csv'),header=True, index=False)
    show_plot = False
    if show_plot : 
        x = np.linspace(0,1)
        y = x**2
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, y, color='blue')
        fig.savefig(os.path.join(output_dir, 'test-plot.png'))
    
    

    html = q2templates.df_to_html(beta, index=False)
    context = {
        'dico': {
            'un': 1, 'deux':2
        },
        'result': html,
        'n_features':len(beta),
        'beta' : beta,
        'show_plot': show_plot,
        'tabs': [{'title': 'Overview',
                  'url': 'overview.html'},
                 {'title': 'LAM fixed',
                  'url': 'lam-fixed.html'}],
        'dangers': [],
        'warnings': [],
    }
    
    
    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    overview_template = os.path.join(TEMPLATES, 'assets', 'overview.html')
    quality_template = os.path.join(TEMPLATES, 'assets', 'quality-plot.html')
    templates = [index, overview_template, quality_template]
    q2templates.render(templates, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'dist'),
                    os.path.join(output_dir, 'dist'))

    


    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("app.init(")
        json.dump({'selected param' : 10}, fh)
        fh.write(',')
        beta.to_json(fh)
        fh.write(');')


'''
Example from : 
https://github.com/qiime2/q2-demux/blob/master/q2_demux/_summarize/_visualizer.py



import collections
import os
import pkg_resources
import shutil
import random
import json

import pandas as pd
import seaborn as sns
import numpy as np

from q2_demux._demux import _read_fastq_seqs
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_demux', '_summarize')


def _decode_qual_to_phred33(qual_str):
    # this function is adapted from scikit-bio
    qual = np.frombuffer(qual_str.encode('ascii'), dtype=np.uint8) - 33
    return qual


# TODO: Remove _PlotQualView once QIIME 2 #220 completed
class _PlotQualView:
    """
    A very simple pass-through view which is made up of a single-end or
    paired-end directory format with a bool indicating if single or paired.
    """
    def __init__(self, directory_format, paired):
        self.directory_format = directory_format
        self.paired = paired


def _link_sample_n_to_file(file_records, counts, subsample_ns):
    results = collections.defaultdict(list)
    for num in subsample_ns:
        total = 0
        for file, sample_id in file_records:
            total += counts[sample_id]
            if num < total:
                idx = counts[sample_id] - (total - num)
                results[file].append(idx)
                break
    return results


def _subsample_paired(fastq_map):
    qual_sample = collections.defaultdict(list)
    min_seq_len = {'forward': float('inf'), 'reverse': float('inf')}
    for fwd, rev, index in fastq_map:
        file_pair = zip(_read_fastq_seqs(fwd), _read_fastq_seqs(rev))
        for i, (fseq, rseq) in enumerate(file_pair):
            if i == index[0]:
                min_seq_len['forward'] = min(min_seq_len['forward'],
                                             len(fseq[1]))
                min_seq_len['reverse'] = min(min_seq_len['reverse'],
                                             len(rseq[1]))
                qual_sample['forward'].append(_decode_qual_to_phred33(fseq[3]))
                qual_sample['reverse'].append(_decode_qual_to_phred33(rseq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    return qual_sample, min_seq_len


def _subsample_single(fastq_map):
    qual_sample = collections.defaultdict(list)
    min_seq_len = {'forward': float('inf'), 'reverse': None}
    for file, index in fastq_map:
        for i, seq in enumerate(_read_fastq_seqs(file)):
            if i == index[0]:
                min_seq_len['forward'] = min(min_seq_len['forward'],
                                             len(seq[1]))
                qual_sample['forward'].append(_decode_qual_to_phred33(seq[3]))
                index.pop(0)
                if len(index) == 0:
                    break
    return qual_sample, min_seq_len


def _compute_stats_of_df(df):
    df_stats = df.describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = df_stats.index.isin(['std', 'mean', 'min', 'max'])
    df_stats = df_stats[~drop_cols]
    return df_stats


def _build_seq_len_table(qscores: pd.DataFrame) -> str:
    sequence_lengths = qscores.notnull().sum(axis=1).copy()
    stats = _compute_stats_of_df(sequence_lengths)

    stats[stats.index != 'count'] = \
        stats[stats.index != 'count'].astype(int).apply('{} nts'.format)

    stats.rename(index={'50%': '50% (Median)',
                        'count': 'Total Sequences Sampled'},
                 inplace=True)
    frame = stats.to_frame(name="")
    return q2templates.df_to_html(frame)


def summarize(output_dir: str, data: _PlotQualView, n: int = 10000) -> None:
    paired = data.paired
    data = data.directory_format
    dangers = []
    warnings = []

    manifest = pd.read_csv(os.path.join(str(data), data.manifest.pathspec),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(data), x))

    fwd = manifest[manifest.direction == 'forward'].filename.tolist()
    rev = manifest[manifest.direction == 'reverse'].filename.tolist()

    per_sample_fastq_counts = {}
    reads = rev if not fwd and rev else fwd
    file_records = []
    for file in reads:
        count = 0
        for seq in _read_fastq_seqs(file):
            count += 1
        sample_id = manifest.loc[manifest.filename == file,
                                 'sample-id'].iloc[0]
        per_sample_fastq_counts[sample_id] = count
        file_records.append((file, sample_id))

    result = pd.Series(per_sample_fastq_counts)
    result.name = 'Sequence count'
    result.index.name = 'Sample name'
    result.sort_values(inplace=True, ascending=False)
    result.to_csv(os.path.join(output_dir, 'per-sample-fastq-counts.csv'),
                  header=True, index=True)
    sequence_count = result.sum()

    if n > sequence_count:
        n = sequence_count
        warnings.append('A subsample value was provided that is greater than '
                        'the amount of sequences across all samples. The plot '
                        'was generated using all available sequences.')

    subsample_ns = sorted(random.sample(range(sequence_count), n))
    link = _link_sample_n_to_file(file_records,
                                  per_sample_fastq_counts,
                                  subsample_ns)
    if paired:
        sample_map = [(file, rev[fwd.index(file)], link[file])
                      for file in link]
        quality_scores, min_seq_len = _subsample_paired(sample_map)
    else:
        sample_map = [(file, link[file]) for file in link]
        quality_scores, min_seq_len = _subsample_single(sample_map)

    forward_scores = pd.DataFrame(quality_scores['forward'])
    forward_stats = _compute_stats_of_df(forward_scores)
    forward_stats.to_csv(os.path.join(output_dir,
                         'forward-seven-number-summaries.csv'),
                         header=True, index=True)
    forward_length_table = _build_seq_len_table(forward_scores)

    if (forward_stats.loc['50%'] > 45).any():
        dangers.append('Some of the PHRED quality values are out of range. '
                       'This is likely because an incorrect PHRED offset '
                       'was chosen on import of your raw data. You can learn '
                       'how to choose your PHRED offset during import in the '
                       'importing tutorial.')

    # Required initilization for conditional display of the table
    reverse_length_table = None
    if paired:
        reverse_scores = pd.DataFrame(quality_scores['reverse'])
        reverse_stats = _compute_stats_of_df(reverse_scores)
        reverse_stats.to_csv(os.path.join(output_dir,
                             'reverse-seven-number-summaries.csv'),
                             header=True, index=True)
        reverse_length_table = _build_seq_len_table(reverse_scores)

    show_plot = len(fwd) > 1
    if show_plot:
        ax = sns.distplot(result, kde=False)
        ax.set_xlabel('Number of sequences')
        ax.set_ylabel('Frequency')
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.png'))
        fig.savefig(os.path.join(output_dir, 'demultiplex-summary.pdf'))

    html_df = result.to_frame().reset_index(drop=False)
    html = q2templates.df_to_html(html_df, index=False)
    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    overview_template = os.path.join(TEMPLATES, 'assets', 'overview.html')
    quality_template = os.path.join(TEMPLATES, 'assets', 'quality-plot.html')
    context = {
        'result_data': {
            'min': result.min(),
            'median': result.median(),
            'mean': result.mean(),
            'max': result.max(),
            'sum': sequence_count
        },
        'forward_length_table': forward_length_table,
        'reverse_length_table': reverse_length_table,
        'result': html,
        'n_samples': result.count(),
        'show_plot': show_plot,
        'paired': paired,
        'tabs': [{'title': 'Overview',
                  'url': 'overview.html'},
                 {'title': 'Interactive Quality Plot',
                  'url': 'quality-plot.html'}],
        'dangers': dangers,
        'warnings': warnings,
    }
    templates = [index, overview_template, quality_template]
    q2templates.render(templates, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'dist'),
                    os.path.join(output_dir, 'dist'))

    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("app.init(")
        json.dump({'n': int(n), 'totalSeqCount': int(sequence_count),
                   'minSeqLen': min_seq_len}, fh)
        fh.write(',')
        forward_stats.to_json(fh)
        if paired:
            fh.write(',')
            reverse_stats.to_json(fh)
        fh.write(');')
'''