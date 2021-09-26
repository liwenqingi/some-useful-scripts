# -*- coding:utf-8 -*-

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

if __name__ == "__main__":
    polyA_cell_lane1_f = "/home/liwenqing/1_work/1_project/1_xiezhe_polyA/lane1_PolyA_result/lane1.polyA.read.count"
    nonpolyA_cell_lane1_f = "/home/liwenqing/1_work/1_project/1_xiezhe_polyA/lane1_PolyA_result/lane1.non_polyA.read.count"
    polyA_cell_lane2_f = "/home/liwenqing/1_work/1_project/1_xiezhe_polyA/lane2_PolyA_result/lane2.polyA.read.count"
    nonpolyA_cell_lane2_f = "/home/liwenqing/1_work/1_project/1_xiezhe_polyA/lane2_PolyA_result/lane2.non_polyA.read.count"
    polyA_4T1RNA_12T_H2_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/polyA_4T1RNA_12T_H2.count"
    nonpolyA_4T1RNA_12T_H2_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/nonpolyA_4T1RNA_12T_H2.count"
    polyA_4T1RNA_12T_SH_H3_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/polyA_4T1RNA_12T_SP_H3.count"
    nonpolyA_4T1RNA_12T_SH_H3_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/nonpolyA_4T1RNA_12T_SH_H3.count"
    polyA_4T1RNA_18T_SP_H3_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/polyA_4T1RNA_18T_SP_H3.count"
    nonpolyA_4T1RNA_18T_SP_H3_f = "/home/liwenqing/2_work/1_project/3_xiezhe_polyA/1_2_Xiezhe_polyA/2_xiezhe_polyA_aligned/nonpolyA_4T1RNA_18T_SP_H3.count"
    
    polyA_cell_lane1_f = pd.read_csv(polyA_cell_lane1_f, sep="\t", index_col=["Geneid"], skiprows=1)
    polyA_cell_lane2_f = pd.read_csv(polyA_cell_lane2_f, sep="\t", index_col=["Geneid"], skiprows=1)
    polyA_4T1RNA_12T_H2_f = pd.read_csv(polyA_4T1RNA_12T_H2_f, sep="\t", index_col=["Geneid"], skiprows=1)
    polyA_4T1RNA_12T_SH_H3_f = pd.read_csv(polyA_4T1RNA_12T_SH_H3_f, sep="\t", index_col=["Geneid"], skiprows=1)
    polyA_4T1RNA_18T_SP_H3_f = pd.read_csv(polyA_4T1RNA_18T_SP_H3_f, sep="\t", index_col=["Geneid"], skiprows=1)
    
    nonpolyA_cell_lane1_f = pd.read_csv(nonpolyA_cell_lane1_f, sep="\t", index_col=["Geneid"], skiprows=1)
    nonpolyA_cell_lane2_f = pd.read_csv(nonpolyA_cell_lane2_f, sep="\t", index_col=["Geneid"], skiprows=1)
    nonpolyA_4T1RNA_12T_H2_f = pd.read_csv(nonpolyA_4T1RNA_12T_H2_f, sep="\t", index_col=["Geneid"], skiprows=1)
    nonpolyA_4T1RNA_12T_SH_H3_f = pd.read_csv(nonpolyA_4T1RNA_12T_SH_H3_f, sep="\t", index_col=["Geneid"], skiprows=1)
    nonpolyA_4T1RNA_18T_SP_H3_f = pd.read_csv(nonpolyA_4T1RNA_18T_SP_H3_f, sep="\t", index_col=["Geneid"], skiprows=1)

    
    polyA_cell_lane1_f = polyA_cell_lane1_f[polyA_cell_lane1_f.columns[-1]]
    polyA_cell_lane2_f = polyA_cell_lane2_f[polyA_cell_lane2_f.columns[-1]]
    polyA_4T1RNA_12T_H2_f = polyA_4T1RNA_12T_H2_f[polyA_4T1RNA_12T_H2_f.columns[-1]]
    polyA_4T1RNA_12T_SH_H3_f = polyA_4T1RNA_12T_SH_H3_f[polyA_4T1RNA_12T_SH_H3_f.columns[-1]]
    polyA_4T1RNA_18T_SP_H3_f = polyA_4T1RNA_18T_SP_H3_f[polyA_4T1RNA_18T_SP_H3_f.columns[-1]]
    
    nonpolyA_cell_lane1_f = nonpolyA_cell_lane1_f[nonpolyA_cell_lane1_f.columns[-1]]
    nonpolyA_cell_lane2_f = nonpolyA_cell_lane2_f[nonpolyA_cell_lane2_f.columns[-1]]
    nonpolyA_4T1RNA_12T_H2_f = nonpolyA_4T1RNA_12T_H2_f[nonpolyA_4T1RNA_12T_H2_f.columns[-1]]
    nonpolyA_4T1RNA_12T_SH_H3_f = nonpolyA_4T1RNA_12T_SH_H3_f[nonpolyA_4T1RNA_12T_SH_H3_f.columns[-1]]
    nonpolyA_4T1RNA_18T_SP_H3_f = nonpolyA_4T1RNA_18T_SP_H3_f[nonpolyA_4T1RNA_18T_SP_H3_f.columns[-1]]
    
    polyA_merge_df = pd.concat([polyA_cell_lane1_f, polyA_cell_lane2_f, polyA_4T1RNA_12T_H2_f, polyA_4T1RNA_12T_SH_H3_f, polyA_4T1RNA_18T_SP_H3_f], axis = 1, join="inner")
    polyA_merge_df.columns = ["polyA_cell_lane1", "polyA_cell_lane2", "polyA_4T1RNA_12T_H2", "polyA_4T1RNA_12T_SP_H3", "polyA_4T1RNA_18T_SP_H3"]
    polyA_merge_df = np.log1p(polyA_merge_df)
    
    nonpolyA_merge_df = pd.concat([nonpolyA_cell_lane1_f, nonpolyA_cell_lane2_f, nonpolyA_4T1RNA_12T_H2_f, nonpolyA_4T1RNA_12T_SH_H3_f, nonpolyA_4T1RNA_18T_SP_H3_f], axis = 1, join="inner")
    nonpolyA_merge_df.columns = ["nonpolyA_cell_lane1", "nonpolyA_cell_lane2", "nonpolyA_4T1RNA_12T_H2", "nonpolyA_4T1RNA_12T_SP_H3", "nonpolyA_4T1RNA_18T_SP_H3"]
    nonpolyA_merge_df = np.log1p(nonpolyA_merge_df)
    
    def corr(x, y, **kwargs):
        coef = stats.spearmanr(x, y)[0]
        label = r"$\rho$ = %.2f"%(coef)
        ax = plt.gca()
        ax.annotate(label, xy=(0.2, 0.9), xycoords=ax.transAxes)
        
    plt.style.use("ggplot")
    g = sns.pairplot(polyA_merge_df, corner=True)
    g.map_lower(corr)
    g.savefig("./pairplot_polyA_correlation.png", dpi=150)
    
    g = sns.pairplot(nonpolyA_merge_df, corner=True)
    g.map_lower(corr)
    g.savefig("./pairplot_nonpolyA_correlation.png", dpi=150)
    
    
    
    
    
    """
    ax1 = fig.add_subplot(1,3,1)
    ax1.scatter(x=merge_lane["lane1"], y = merge_lane["lane2"])
    rho, p = stats.spearmanr(merge_lane["lane1"], merge_lane["lane2"])
    ax1.text(s=r"$\rho$: %.2f"%(rho),x=0.7, y=0.7)
    ax1.set_xlabel("Lane1 log1p(gene count)")
    ax1.set_ylabel("Lane2 log1p(gene count)")
    
    ax2 = fig.add_subplot(1,3,2)
    ax2.scatter(x=merge_lane["lane1"], y = merge_lane["lane3"])
    rho, p = stats.spearmanr(merge_lane["lane1"], merge_lane["lane3"])
    ax2.text(s=r"$\rho$: %.2f"%(rho),x=0.7, y=0.7)
    ax2.set_xlabel("Lane1 log1p(gene count)")
    ax2.set_ylabel("Lane3 log1p(gene count)")
    
    ax3 = fig.add_subplot(1,3,3)
    ax3.scatter(x=merge_lane["lane2"], y = merge_lane["lane3"])
    rho, p = stats.spearmanr(merge_lane["lane2"], merge_lane["lane3"])
    ax3.text(s=r"$\rho$: %.2f"%(rho),x=0.7, y=0.7)
    ax3.set_xlabel("Lane2 log1p(gene count)")
    ax3.set_ylabel("Lane3 log1p(gene count)")
    
    fig.savefig("./scatter_correlation_nonpolyA.png",dpi=100)
    """