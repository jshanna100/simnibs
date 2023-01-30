import pickle
import pandas as pd
from os.path import join
import seaborn as sns
import matplotlib.pyplot as plt
plt.ion()

root_dir = "/home/jev/"

filepath_3 = join(root_dir, "simnibs", "3_emp", "success_record.pickle")
filepath_4 = join(root_dir, "simnibs", "4_emp", "success_record.pickle")

with open(filepath_3, "rb") as f:
    res_3 = pickle.load(f)
with open(filepath_4, "rb") as f:
    res_4 = pickle.load(f)

dfs = {}
for res_name, res in zip(["3", "4"], [res_3, res_4]):
    df_dict = {"Subj":[], "Proj":[], "Result":[]}
    for r in res:
        df_dict["Subj"].append(r[0])
        df_dict["Proj"].append(r[1])
        df_dict["Result"].append(r[2])
    dfs[res_name] = pd.DataFrame.from_dict(df_dict)
