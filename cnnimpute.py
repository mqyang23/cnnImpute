import rpy2.robjects as robjects
import concurrent.futures
from tensorflow.python.client import device_lib
import timeit
from sklearn.preprocessing import Normalizer
import numpy as np
import pandas as pd
from cnnImpute.multimodel import *
from cnnImpute.parser import parse_args
from functools import partial
#inorder to pass args to r, import this
import subprocess

def get_available_gpus():
  local_device_protos = device_lib.list_local_devices()
  return [x.name for x in local_device_protos if x.device_type == 'GPU']



def cnnImpute(**kwargs):
  args = parse_args()
  for key, value in kwargs.items():
    setattr(args, key, value)
  start_begin = timeit.default_timer()
  print("Start labelling data...")
  r_script = "cnnImpute/R_code/clusteringdata_missingvalue.r"
  command = ['Rscript', r_script, args.inputFile]
  subprocess.call(command, stderr=subprocess.STDOUT)
  normalization = False
  tmpfile = args.inputFile.rsplit('/',1)[1]
  data = pd.read_csv("cnnImpute/tmp_file/" + tmpfile, index_col=0)
  print("Overview of labeled data:")
  print(data.iloc[0:5,0:5])
  print("data preprocess...")
  geneall0 = []
  count = np.array(data.sum(axis=1))

  for row in range(0, data.shape[0]):
    if count[row] == 0:
      geneall0.append(row)

  index = np.array(data.index)
  data = data.T
  # convert data frame to array
  index1 = np.array(data.index)

  data_original = data.copy()
  data_original = data_original.astype(float)
  data = np.log1p(data)

  # log normalize data
  data = np.array(data)
  data1 = data.copy()

  if normalization == True:
    data_rescaled = pd.DataFrame(Normalizer().fit_transform(data_original))
    data_rescaled = np.log1p(data_rescaled)
  else:
    data_rescaled = pd.DataFrame(data)

  # extract cluster information into index1
  for i in range(0, len(index1)):
    # index1[i]=str(index1[i]).split("_")[0].split("'")[1]
    index1[i] = str(index1[i]).split("_")[0]
  index1 = index1.astype(int)

  # seperate index into different cluster
  for i in range(1, int(max(index1)) + 1):
    text = "cluster" + str(i)  # str means string
    globals()[str(text)] = np.where(index1 == i)[0]

  # save file into csv
  to_csvfile = data_original.transpose()
  print("write file:", to_csvfile.iloc[0:4, 0:4])

  to_csvfile.to_csv("cnnImpute/tmp_file/data_original.csv")

  robjects.r('source("cnnImpute/R_code/scimpute_v3.R")')

  #################################################################################################################333
  clust = []

  gpus = get_available_gpus()
  gpu_Cnt = len(gpus)
  global ind
  k = max(index1)
  with concurrent.futures.ThreadPoolExecutor(len(gpus)) as executor:
    for ind in range(0, max(index1), gpu_Cnt):
      if k > 0:
        partial_res = partial(runModel, args.learning_rate, args.batch_size,
                                   args.max_epochs, args.dropout_rate, args.dropoutRate_threshold,
                                   args.inputGene_threshold,data,data1,data_rescaled,ind,geneall0,index1)
        if k >= gpu_Cnt:
          results = [x for x in executor.map(partial_res, gpus)]
        else:
          gpus = gpus[0:k]
          results = [x for x in executor.map(partial_res, gpus)]
        results = pd.concat(results, axis=0, ignore_index=True)
        if ind == 0:
          outputs = results
        else:
          outputs = pd.concat([outputs, results], axis=0)
        k -= gpu_Cnt
  for i in range(1, int(max(index1)) + 1):
      text = "cluster" + str(i)  # str means string
      globals()[str(text)] = np.where(index1 == i)[0]

  for i in range(max(index1)):
    clust += list(globals()["cluster" + str(i + 1)])

  res = outputs.copy()

  for rowx in range(0, res.shape[0]):
    loc = clust[rowx]
    res.iloc[loc, :] = outputs.iloc[rowx, :]
  res = np.expm1(res)
  res = pd.DataFrame(res)
  print("output file： ",args.output)
  res.to_csv(args.output)
  stop_end = timeit.default_timer()
  print('total Time: ', stop_end - start_begin)