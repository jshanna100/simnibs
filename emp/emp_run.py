from simnibs import __version__
from simu_func import emp_montage
from itertools import product

version = int(__version__[0])

projs = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
#root_dir = "/home/jev/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))
subj_dicts = build_subject_paths(data_dir, version)
n_jobs = 6

queue = list(product(subj_dicts, projs))
Parallel(n_jobs=n_jobs)(delayed(emp_montage)(*q, root_dir) for q in queue)
