import sys

import numpy as np
import pandas as pd
import plotly.express as px

def load_tsv(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t', header=0, index_col=False)
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