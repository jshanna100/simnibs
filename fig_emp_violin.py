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
root_dir = "/home/jev/"

fig_dir = join(root_dir, "simnibs/figures")

# data_dir = join(root_dir, "simnibs/3_emp")
# df_3 = pd.read_pickle(join(data_dir, "df_emp_3.pickle"))
# df_3["Version"] = ["3"] * len(df_3)
data_dir = join(root_dir, "simnibs/4_emp")
df_4 = pd.read_pickle(join(data_dir, "df_emp_4.pickle"))
df_4["Version"] = ["4"] * len(df_4)
#df = pd.concat([df_3, df_4])


sns.violinplot(data=df_4, x="Project", y="Mag", hue="Version",
               order=[f"P{p}" for p in np.arange(1,9)], inner="points")
