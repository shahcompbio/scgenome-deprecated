import pandas as pd
import numpy as np
from scgenome import qc
from os.path import join
import os
import sys

arg = {
  "out_dir": "/work/shah/maherm/tantalus/SC-2534",
  "cn_data_fp": "/work/shah/maherm/tantalus/SC-2534/A96200B_reads.csv",
  "metrics_fp": "/work/shah/maherm/tantalus/SC-2534/metrics.csv",
}

print(f"sys.argv {sys.argv}")
if len(sys.argv) >= 3:
    arg["out_dir"] = sys.argv[1]
    arg["cn_data_fp"] = sys.argv[2]
    arg["metrics_fp"] = sys.argv[3]
if not os.path.exists(arg["out_dir"]):
    print(arg["out_dir"] + " does not exist, creating it")
    os.makedirs(arg["out_dir"])
pd.DataFrame(arg, index=[0]).to_csv(join(arg["out_dir"], "qc_arg.csv"), 
                                    index=None)

cn_data = pd.read_csv(arg["cn_data_fp"])
metrics = pd.read_csv(arg['metrics_fp'])

qcn, qcn_data = qc.qc_cn(metrics, cn_data, quality=None)
qcn_data.to_csv(join(arg["out_dir"], "qc_cn_data.csv"), index=None)


