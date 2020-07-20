import pandas as pd
import os
import random
from collections import defaultdict
from subprocess import PIPE, Popen


def merge(df):
    results = []
    cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = \
        None, None, None, None, None, None, None, None
    for i in range(df.shape[0]):
        chr, start, end, center, width_above_cutoff, total_signal, height, height_logP = df.ix[i, :]
        if cur_chr is None:
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        elif cur_chr != chr:
            results.append([cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP])
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        elif start - cur_end > 3000:
            results.append(
                [cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height,
                 cur_height_logP])
            cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height, cur_height_logP = chr, start, end, center, width_above_cutoff, total_signal, height, height_logP
        else:
            cur_end = end
            cur_center = (cur_start + cur_end)/2
            cur_width_above_cutoff += width_above_cutoff
            cur_total_signal += total_signal
            cur_height = cur_height if height < cur_height else height
            cur_height_logP = cur_height_logP if height_logP < cur_height_logP else height_logP
    results.append([cur_chr, cur_start, cur_end, cur_center, cur_width_above_cutoff, cur_total_signal, cur_height,
                 cur_height_logP])
    final_df = pd.DataFrame(results)
    final_df.columns = df.columns
    return final_df


def random_peaks(peaks, size, seed):
    """

    :param peaks: list of peaks with start, end, chr
    :param size: chromosome size, determine the size of random array
    :param seed: random seed
    :return: a random peaks with similar size
    """
    random.seed(seed)

    max_width = max([peak[4] for peak in peaks])

    new_starts = random.sample(range(1, size-max_width), len(peaks))
    new_starts = sorted(new_starts)
    peaks = sorted(peaks, key=lambda x:x[4], reverse=True)

    new_peaks = []

    for i in range(len(peaks)):
        chr, start, end, center, width_above_cutoff, total_signal, height, height_logP = peaks[i]
        new_start = new_starts[i]
        new_end = new_start + (end - start)
        new_peaks.append([chr,
                          new_start, new_end,
                          center, width_above_cutoff,
                          total_signal, height, height_logP])
    return new_peaks


def sort_bed(danpos_result):
    df = pd.read_csv(danpos_result, sep='\t')
    df.to_csv(danpos_result.replace(".xls", ".bed"), sep='\t', index=None, header=None)

    command = ["sort", "-k1,1", "-k2,2n", danpos_result, ">", danpos_result.replace('.xls', '_sorted.bed')]
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    _, error = process.communicate()

    return danpos_result.replace('.xls', '_sorted.bed')


def merge_bed(bed1, bed2, out_name="merged_feature.bed"):
    command = ["bedtools", "closest", "-a", bed1, "-b", bed2, ">", out_name]
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    _, error = process.communicate()

    overlap_df = pd.read_csv(out_name, sep='\t', header=None)

    results = []
    for k in range(overlap_df.shape[0]):
        chr, start1, end1, start2, end2 = overlap_df.iloc[k, [0, 1, 2, 9, 10]]
        if start1 <= start2 <= end1 or start1 <= end2 <= end1 or start2 <= start1 <= end2 or start2 <= end1 <= end2:
            distance = 0
        else:
            distance = abs(end1 - start2) if abs(end1 - start2) < abs(start1 - end2) else abs(start1 - end2)
        results.append([chr, start1, end1, start2, end2, distance])
    final_df = pd.DataFrame(results)
    final_df.columns = ['chr', 'start1', 'end1',
                        'start2', 'end2',
                        'distance']
    return final_df


# merge GTF and combined two features bed, and calculate the number of overlap
def calculate_overlap(gene_table, gene_bed, bed_peak):
    merged_df = merge_bed(gene_bed, bed_peak)
    merged_df = merged_df[merged_df['distance'] == 0]
    gene_bed_df = pd.read_csv(gene_bed, sep='\t')
    merged_df = merged_df.set_index(["chr", "start", "end"])
    gene_bed_df = gene_bed_df.set_index(["chr", "start", "end"])

    gene_bed_df = gene_bed_df[gene_bed_df.index.isin(merged_df.index)]
    genes = gene_bed_df["hg19.kgXref.geneSymbol"].unique()

    gene_table.loc[genes, "overlap"] = 1

    return gene_table