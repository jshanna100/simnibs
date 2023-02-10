import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.gridspec as gridspec
import seaborn as sns
import os
import numpy as np
from os.path import join
from scipy.stats import mode

import matplotlib
font = {'weight' : 'bold',
        'size'   : 18}
matplotlib.rc('font', **font)

root_dir = "/home/hannaj/"
#root_dir = "/home/jev/"

fig_dir = join(root_dir, "simnibs/figures")

data_dir = join(root_dir, "simnibs/3_emp")
df_3 = pd.read_pickle(join(data_dir, "df_emp_3.pickle"))
df_3["Version"] = ["3"] * len(df_3)
data_dir = join(root_dir, "simnibs/4_emp")
df_4 = pd.read_pickle(join(data_dir, "df_emp_4.pickle"))
df_4["Version"] = ["4"] * len(df_4)
df = pd.concat([df_3, df_4])

df = df.query("Summary=='ROI_Median'")
order = ["P1", "P2",  "P3", "P4", "P5", "P7", "P8"]
facet = sns.catplot(data=df, x="Project", y="Mag", hue="Version", kind="violin",
                    order=order, inner=None)
facet = sns.stripplot(data=df, x="Project", y="Mag", hue="Version",
                    order=order, dodge=True, color="black", legend=False)
facet.axes.axhline(df.query("Version=='3' and Summary=='ROI_Median'")["Mag"].mean(),
                    color="tab:blue", linestyle="--")
facet.axes.axhline(df.query("Version=='4' and Summary=='ROI_Median'")["Mag"].mean(),
                    color="tab:orange", linestyle="--")

# facet = sns.catplot(data=df, x="Project", y="Foc", hue="Version", kind="violin",
#                     order=order, inner="points")
facet = sns.catplot(data=df, x="Project", y="Foc", hue="Version", kind="violin",
                    order=order, inner=None)
facet = sns.stripplot(data=df, x="Project", y="Foc", hue="Version",
                    order=order, dodge=True, color="black", legend=False)
