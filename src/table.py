import os
import pandas as pd
from collections import defaultdict
from multiprocessing import Process, Queue
from util import df_to_index_sk, df_to_index_danpos


def get_range_absolute(gene_list, all_gene_GTF, left_distance, right_distance, TSS_pos, TTS_pos):
    """
    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """

    if gene_list is not None:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'].isin(gene_list)]
    else:
        cur_df = all_gene_GTF

    positive_df = cur_df[cur_df['hg19.knownGene.strand'] == '+'].copy()
    negative_df = cur_df[cur_df['hg19.knownGene.strand'] == '-'].copy()

    if TSS_pos == 'TSS' and TTS_pos == 'TSS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txStart'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txEnd'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'TSS' and TTS_pos == 'TTS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txEnd'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txStart'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'MID':
        positive_df['MID'] = (positive_df['hg19.knownGene.txStart'] + positive_df['hg19.knownGene.txEnd'])/2
        negative_df['MID'] = (negative_df['hg19.knownGene.txStart'] + negative_df['hg19.knownGene.txEnd'])/2

        positive_df['left_range'] = positive_df['MID'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['MID'] + right_distance

        negative_df['right_range'] = negative_df['MID'] - left_distance
        negative_df['left_range'] = negative_df['MID'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0

    positive_df['length'] = positive_df['hg19.knownGene.txEnd'] - positive_df['hg19.knownGene.txStart']
    negative_df['length'] = negative_df['hg19.knownGene.txEnd'] - negative_df['hg19.knownGene.txStart']

    new_df = positive_df.append(negative_df)

    if len(new_df['hg19.kgXref.geneSymbol'].unique()) != len(all_gene_GTF['hg19.kgXref.geneSymbol'].unique()):
        new_df.to_csv('wrong_range.csv')

    result_df = new_df[['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'left_range', 'right_range', 'length']]
    result_df.columns = ['gene', 'chr', 'left_range', 'right_range', 'length']

    return result_df


def get_stats(gene_df, df_path, criteria, cell_type, bin=3000, df_function=df_to_index_danpos):
    """
    This function will select the target gene's peaks stats from the dataframe
    :param gene_df: gene:range dictionary, get the best result for each gene from transcripts range
    :param df: dataframe contain peaks from certain cutoff
    :return:
    """
    # print 'get stats', len(gene_df['gene'].unique())
    if criteria != 'skewness' and criteria != 'kurtosis':
        table_dict = df_function(df_path)
    else:
        df_function = df_to_index_sk
        table_dict = df_function(df_path)

    results = defaultdict(float)

    for k in range(gene_df.shape[0]):
        gene_name = gene_df.iloc[k, 0]

        chr_name, start, end, length = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3], gene_df.iloc[k, 4]
        ## Here is the problem, danpos selector will consider the entire overlapped peaks
        ## The other approach is using self designed peak calling, to make sure each parameter will return different value
        cur_table = set()

        if end < start:
            mid = (start + end) / 2
            start = mid
            end = mid

        for i in range(int(start/bin), int(end/bin) + 1):
            if chr_name in table_dict and i in table_dict[chr_name]:
                table = table_dict[chr_name][i]
                cur_table = cur_table.union(table)

        if len(cur_table) == 0:
            continue

        selected_table = []
        for t in cur_table:
            if start < t[1] < end:
                selected_table.append(t)
            elif start < t[2] < end:
                selected_table.append(t)
            elif t[1] <= start and end <= t[2]:
                selected_table.append(t)

        if len(selected_table) == 0:
            continue

        cur_df = pd.DataFrame(list(selected_table))

        if cur_df.shape[1] == 6:
            cur_df.columns = ['chr',
                          'start',
                          'end',
                          'width_above_cutoff',
                          'total_signal',
                          'height',]
        else:
            cur_df.columns = ['chr',
                              'start',
                              'end',
                              'width_above_cutoff',
                              'total_signal',
                              'height',
                              'skewness',
                              'kurtosis']

        if criteria == 'total_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.sum()
        elif criteria == 'height':
            cur_value = cur_df['height'].max()
        elif criteria == 'single_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.max()
        elif criteria == 'total_signal':
            cur_value = cur_df['total_signal'].sum()
        elif criteria == 'single_signal':
            cur_value = cur_df['total_signal'].max()
        elif criteria == 'coverage':
            cur_value = (cur_df['end'] - cur_df['start']).sum()*1.0/length

        # # This is for kurtosis and skewness
        elif cur_df.shape[0] > 0 and criteria == 'skewness' and 'skewness' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(),'skewness']
        elif cur_df.shape[0] > 0 and criteria == 'kurtosis' and 'kurtosis' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(), 'kurtosis']


        if cur_value > results[gene_name] and criteria != 'skewness' and criteria != 'kurtosis':
            results[gene_name] = cur_value
        # this is for kurtosis and skewness

        elif criteria == 'kurtosis':
            if abs(cur_value) > abs(results[gene_name]):
                results[gene_name] = cur_value
        elif criteria == 'skewness':
            if abs(cur_value) > results[gene_name]:
                results[gene_name] = abs(cur_value)

    final = []

    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name], cell_type))
    print len(final)
    return final


def selector(all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria, cell_type,
                     TSS_pos, TTS_pos):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    all_gene_ranges = get_range_absolute(None, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos, TTS_pos)

    cur_df = all_dfs[cell_type][cutoff]

    # overlap selecter
    all_gene_results = get_stats(all_gene_ranges, cur_df, criteria, cell_type)

    all_gene_results_df = pd.DataFrame(all_gene_results)

    all_gene_results_df.columns = ['gene', criteria, 'cell_type']
    all_gene_results_df = all_gene_results_df[['gene', criteria]]
    return all_gene_results_df


ref_data_root = "../reference/"

def get_real_table(target_cell_type, peaks_path, out_name, gtf_path=ref_data_root+"hg19GREATgene2UCSCknownGenes.table.xls",
                   best_parameter_path=ref_data_root+'parameters_default.csv',
                   process=8):
    ## get the CIG genesets and nonCIG genesets
    all_gene_GTF = pd.read_csv(gtf_path, sep='\t', dtype={'hg19.kgXref.geneSymbol': str})

    # all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    # get feature parameter table
    best_feature_df = pd.read_csv(best_parameter_path)

    # get all tables, if want to use different parameters for different markers
    all_tables = {}

    # all peak_files from same cell type will be in one folder peak_path_root/celltype/marker/
    m_types = os.listdir(peaks_path+'/'+target_cell_type)
    for m in m_types:
        if m not in all_tables:
            all_tables[m] = defaultdict(dict)
            dfs_path = peaks_path+'/'+target_cell_type+'/'+m+'/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]

            for table_name in dfs:
                info = table_name.split('_')
                cell_type = target_cell_type
                cutoff = float(info[-1][:-4])
                all_tables[m][cell_type][cutoff] = dfs_path+table_name

    chunks = []
    cur_index = 0
    reminder = best_feature_df.shape[0] % process
    chunk_size = best_feature_df.shape[0] / process
    for i in range(process):
        if reminder > 0:
            chunks.append(best_feature_df[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
            cur_index += 1
            reminder -= 1
        else:
            chunks.append(best_feature_df[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

    total_chunk_size = 0
    for chunk in chunks:
        total_chunk_size += len(chunk)
    if total_chunk_size != best_feature_df.shape[0]:
        print 'multiple processes chunk size is not correct'
        return None

    queue = Queue()
    processes = []

    for i in range(process):
        cur_chunk = chunks[i]
        p = Process(target=real_table_process,
                    args=(queue, target_cell_type, cur_chunk, all_gene_GTF, all_tables, all_sk_tables, m_types))
        processes.append(p)
        p.start()

    final_feature_tables = []

    for i in range(process):
        cur_final_feature_df = queue.get()
        final_feature_tables.append(cur_final_feature_df)
    final_feature_df = pd.concat(final_feature_tables, axis=1)

    for p in processes:
        p.join()

    index_cols = []
    index_col = False
    for ci in range(len(final_feature_df.columns)):
        if final_feature_df.columns[ci].startswith('gene'):
            if not index_col:
                index_cols.append(ci)
                index_col = True
        else:
            index_cols.append(ci)
    final_feature_df = final_feature_df.iloc[:, index_cols]
    columns = list(final_feature_df.columns)
    columns[0] = 'gene_id'
    final_feature_df.columns = columns

    if out_name is None:
        final_feature_df.to_csv(target_cell_type+'_' + best_parameter_path[best_parameter_path.rfind('/')+1:-4]+'_'+'real_table.csv', index=None)
    else:
        final_feature_df.to_csv(out_name, index=None)
    return final_feature_df


def real_table_process(queue, target_cell_type, best_feature_df, all_gene_GTF, all_tables, all_sk_tables, m_types):
    cur_final_feature_df = None
    # print best_feature_df.columns
    for i in range(best_feature_df.shape[0]):
        marker = best_feature_df.iloc[i, 0]
        if marker not in m_types:
            continue
        feature = best_feature_df.iloc[i, 1]

        criteria = feature.replace('_genebody', '')

        start = best_feature_df.iloc[i, 2]
        end = best_feature_df.iloc[i, 3]
        height = best_feature_df.iloc[i, 4]
        # print marker, feature, start, end, height
        if feature.find('genebody') != -1:
            TSS, TTS = 'TSS', 'TTS'
        else:
            TSS, TTS = 'TSS', 'TSS'

        if feature.find('kurtosis') != -1 or feature.find('skewness') != -1:
            option = True
        else:
            option = False

        if option:
            cur_stat_dfs = all_sk_tables[marker]
        else:
            cur_stat_dfs = all_tables[marker]

        all_gene_results_df = selector(all_gene_GTF, start, end, cur_stat_dfs, height,
                                               criteria, target_cell_type, TSS, TTS)

        all_gene_results_df.columns = ['gene'+str(i), marker + '_' + feature]
        if cur_final_feature_df is None:
            cur_final_feature_df = all_gene_results_df.copy()
        else:
            cur_final_feature_df[marker + '_' + feature] = all_gene_results_df[all_gene_results_df.columns[1]]

    queue.put(cur_final_feature_df)
    return
