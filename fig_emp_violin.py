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
        'size'   : 28}
matplotlib.rc('font', **font)
#sns.set(rc={"figure.figsize":(19.2, 19.2)})


root_dir = "/home/hannaj/"
root_dir = "/home/jev/"

fig_dir = join(root_dir, "simnibs/figures/")

data_dir = join(root_dir, "simnibs/3_emp")
df_3 = pd.read_pickle(join(data_dir, "df_emp_3.pickle"))
df_3["Version"] = ["3"] * len(df_3)
data_dir = join(root_dir, "simnibs/4_emp")
df_4 = pd.read_pickle(join(data_dir, "df_emp_4.pickle"))
df_4["Version"] = ["4"] * len(df_4)
df = pd.concat([df_3, df_4], ignore_index=True)

df = df.query("Summary=='ROI_Median'")
order = ["P1", "P2", "P3", "P4", "P5", "P7", "P8"]

#df = df.query("Version=='3'")
hue_spec = "Version"

dot_size = 9
mag_fig, mag_ax = plt.subplots(1, 1, figsize=(25.6, 14.4))
facet = sns.violinplot(data=df, x="Project", y="Mag", hue=hue_spec,
                    order=order, inner=None, ax=mag_ax)
facet = sns.stripplot(data=df, x="Project", y="Mag", hue=hue_spec,
                      order=order, dodge=True, color="black", legend=False,
                      size=dot_size, ax=mag_ax)
facet.axes.axhline(df.query("Version=='3' and Summary=='ROI_Median'")["Mag"].mean(),
                    color="tab:blue", linestyle="--")
facet.axes.axhline(df.query("Version=='4' and Summary=='ROI_Median'")["Mag"].mean(),
                    color="tab:orange", linestyle="--")
plt.title("Montage magnitudes", fontsize=48)
mag_fig.savefig(f"{fig_dir}montage_mag.pdf")

foc_fig, foc_ax = plt.subplots(1, 1, figsize=(25.6, 14.4))
facet = sns.violinplot(data=df, x="Project", y="Foc", hue=hue_spec,
                    order=order, inner=None, ax=foc_ax)
facet = sns.stripplot(data=df, x="Project", y="Foc", hue=hue_spec,
                    order=order, dodge=True, color="black", legend=False,
                    size=dot_size, ax=foc_ax)
plt.title("Montage focalities", fontsize=48)
foc_fig.savefig(f"{fig_dir}montage_foc.pdf")
