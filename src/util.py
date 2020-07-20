import pandas as pd
import numpy as np


def df_to_index_danpos(df, bin=3000):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        line = line.split()
        if len(line) <=1:
            continue
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results


def df_to_index_sk(df, bin=3000):
    results = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.strip().split(',')
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             float(line[8]),
             float(line[9]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results


def get_header(go_term_table):
    n = 0
    with open(go_term_table, 'r') as f:
        for line in f.readlines():
            if (line.find("GO biological process complete") != -1):
                return n
            n += 1
    return n


def pathway_enrich(pathways, go_term_results):
    n = get_header(go_term_results)

    df = pd.read_csv(f, sep='\t', index_col=0, header=n)

    df = df.loc[pathways, u'upload_1 (raw P-value)']
    return df


def calc_MI(X, Y, bins):
   c_XY = np.histogram2d(X,Y,bins)[0]
   c_X = np.histogram(X,bins)[0]
   c_Y = np.histogram(Y,bins)[0]

   H_X = shan_entropy(c_X)
   H_Y = shan_entropy(c_Y)
   H_XY = shan_entropy(c_XY)

   MI = H_X + H_Y - H_XY
   return MI


def shan_entropy(c):
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log2(c_normalized))
    return H


def MI_matrix(df, features):
    A = df[features].values
    bins = int(df.max().max()/1000)
    n = df.shape[1]
    matMI = np.zeros((n, n))

    for ix in np.arange(n):
        matMI[ix, ix] = np.nan
        for jx in np.arange(ix+1,n):
            matMI[ix, jx] = calc_MI(A[:,ix], A[:,jx], bins)
            matMI[jx, ix] = matMI[ix, jx]

    final_df = pd.DataFrame(matMI)
    final_df.columns = features
    final_df.index = features
    return final_df