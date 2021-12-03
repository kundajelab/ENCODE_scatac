import sys

import numpy as np
import pandas as pd
import plotly.express as px

def load_tsv(tsv_path):
    with open(tsv_path) as f:
        head = [next(f) for _ in range(4)]
        header_data = []
        for l in head:
            curr = None
            entries = []
            for i in l.rstrip("\n").split("\t"):
                if i == "":
                    entries.append(curr)
                else:
                    entries.append(i)
                    curr = i
            header_data.append(entries)
        header = [":".join(i) for i in zip(*header_data)]
        df = pd.read_csv(f, sep='\t', names=header, index_col=False)

    cols_to_use = [i for i, v in df.dtypes.items() if np.issubdtype(v, np.number)]
    return df, cols_to_use

def plot_hist(df, col):
    fig = px.histogram(df, x=col, marginal="rug", hover_data=df.columns)
    return fig

def plot_qc(tsv_path, out_path):
    df, cols = load_tsv(tsv_path)
    with open(out_path, 'a') as f:
        for col in cols:
            h = plot_hist(df, col)
            f.write(h.to_html(full_html=False, include_plotlyjs='cdn'))
        
if __name__ == '__main__':
    tsv_path = sys.argv[1]
    out_path = sys.argv[2]
    plot_qc(tsv_path, out_path)