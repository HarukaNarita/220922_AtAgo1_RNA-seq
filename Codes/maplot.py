from operator import index
from textwrap import indent
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class MAPlot:
    f_fC = "/Users/siomi19/Desktop/siomi_lab/Files/testis_RNA-seq_featureCounts.txt"
    dir_fC_tr = "/Users/siomi19/Desktop/siomi_lab/Files/testis_RNA-seq_featureCounts_triplicate"
    age_to_column = {"E135":0, "E165":3, "P0":6, "P3":9}
    drop = [0,1,2,3,4]
    dir_MA = "/Users/siomi19/Desktop/siomi_lab/Figures"
    dir_DEG = "/Users/siomi19/Desktop/siomi_lab/Files"
    prefix = "AtAgo1_RNA-seq_featureCounts"
    replicates = 3
    gene_type = ""

    def __init__(self, age_control:str, age_treat:str):
        self.age_control = age_control
        self.age_treat = age_treat
        self._df_DESeq2 = None
    
    @property
    def df_DESeq2(self):
        try:
            self._df_DESeq2 = pd.read_csv(f"{MAPlot.dir_fC_tr}/{MAPlot.prefix}_{MAPlot.gene_type}_{self.age_control}_{self.age_treat}_DESeq2.txt",sep="\t", skiprows=0, header=0, index_col=0)
            return self._df_DESeq2
        except FileNotFoundError as e:
            return None
    @df_DESeq2.setter
    def df_DESeq2(self, df):
        self.df_DESeq2 = df

    def process_df(self):
        df_fC = pd.read_csv(MAPlot.f_fC, sep="\t", skiprows=0, header=1,index_col=0)
        df_fC = df_fC.drop(df_fC.columns[MAPlot.drop], axis=1)    # remove outliers to make triplecates

        col_c = MAPlot.age_to_column[self.age_control]
        col_t = MAPlot.age_to_column[self.age_treat]
        
        df_fC_control = df_fC.iloc[:, col_c:col_c + MAPlot.replicates]
        df_fC_treat = df_fC.iloc[:, col_t:col_t + MAPlot.replicates]

        df_fC_control = df_fC_control.set_axis([f"{self.age_control}-{r}" for r in range(1,MAPlot.replicates+1,1)], axis=1)
        df_fC_treat = df_fC_treat.set_axis([f"{self.age_treat}-{r}" for r in range(1,MAPlot.replicates+1,1)], axis=1)
        
        self._df_fC = pd.concat([df_fC_control,df_fC_treat],axis=1)
        self._df_fC.to_csv(f"{MAPlot.dir_fC_tr}/{MAPlot.prefix}_{MAPlot.gene_type}_{self.age_control}_{self.age_treat}.txt", sep="\t", header=True, index=True)
    
    def plot_MA(self, upregulated:bool, q_values:list, fc_threshold:float, bh=True, save=True):
        df = self.df_DESeq2
        plt.figure()
        cmap = plt.get_cmap("rainbow")
        plt.scatter(np.log(df["baseMean"])/np.log(2), df["log2FC"], c="gray", s=1)
        if bh:
            method = "q.value_BH"
        else:
            method = "q.value"
        
        if upregulated:
            lfc = (df["log2FC"] > fc_threshold)
        else:
            lfc = (df["log2FC"] < fc_threshold)

        for i, q_value in enumerate(q_values):
            c = np.array(cmap(i/len(q_values))).reshape(1,-1)
            plt.scatter(np.log(df["baseMean"][df[method]<q_value][lfc])/np.log(2),df["log2FC"][df[method]<q_value][lfc], c=c, s=1, label=f"q value = {q_value}")
        plt.axhline(y=0, c="black")
        plt.xlabel("mean normalized counts")
        plt.ylabel("log fold change")
        plt.title(f"{MAPlot.gene_type} : {self.age_control} vs {self.age_treat} ( log fold change < {fc_threshold} )")
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        if save:
            q_values_name = [str(q) for q in q_values]
            ifBH = lambda : "_BH" if bh else ""
            plt.savefig(f"{MAPlot.dir_MA}/{MAPlot.gene_type}/MA-plot_{self.age_control}_{self.age_treat}_{'_'.join(q_values_name)}_{fc_threshold}{ifBH()}.jpg", bbox_inches="tight")

    def get_DEG(self, upregulated:bool, q_values:list, fc_threshold:float, bh=True, save=True):
        df = self.df_DESeq2
        DEGs = []
        if bh:
            method = "q.value_BH"
        else:
            method = "q.value"

        if upregulated:
            lfc = (df["log2FC"] > fc_threshold)
        else:
            lfc = (df["log2FC"] < fc_threshold)

        for q_value in q_values:
            df_DEG = df[df[method]<q_value][lfc]
            DEGs.append(df_DEG)
        if save:
            ifBH = lambda : "_BH" if bh else ""
            for i, df in enumerate(DEGs):
                df.to_csv(f"{MAPlot.dir_DEG}/{MAPlot.gene_type}/DEG_{self.age_control}_{self.age_treat}_{q_values[i]}_{fc_threshold}{ifBH()}.csv", sep="\t", index=True)
        return DEGs


