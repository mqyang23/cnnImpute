import rpy2.robjects as robjects
import concurrent.futures
from tensorflow.python.client import device_lib
import timeit
from sklearn.preprocessing import Normalizer
import numpy as np
import pandas as pd
from cnnimpute.multimodel import *
from cnnimpute.parser import parse_args
from functools import partial

def get_available_gpus():
  local_device_protos = device_lib.list_local_devices()
  return [x.name for x in local_device_protos if x.device_type == 'GPU']



def cnnImpute(**kwargs):
  args = parse_args()

  for key, value in kwargs.items():
    setattr(args, key, value)
  start_begin = timeit.default_timer()
  print("Start labelling data...")
  robjects.r('source("cnnimpute/R_code/clusteringdata_missingvalue.r")')

  # for thre in [0,0.01,0.1]:
  # for thre in [0]:
  normalization = False
  # address = '/home/junie/recovereddata/'
  # address = '/home/junie/'
  # path = address + filename + ".csv"

  data = pd.read_csv("cnnimpute/tmp_file/"+str(args.inputFile), index_col=0)
  print("Overview of labeled data:")
  print(data.iloc[0:5,0:5])
  print("data preprocess...")
  geneall0 = []
  count = np.array(data.sum(axis=1))

  for row in range(0, data.shape[0]):
    if count[row] == 0:
      geneall0.append(row)

  index = np.array(data.index)

  # transpose data
  # data = data.transpose()
  data = data.T
  # convert data frame to array
  index1 = np.array(data.index)

  # data = data.astype(int)
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
  # print("write file:", to_csvfile.iloc[0:4, 0:4])

  to_csvfile.to_csv("cnnimpute/tmp_file/data_original.csv")
  print("Estimate missing value probability...")
  robjects.r('source("cnnimpute/R_code/scimpute_v3.R")')

  #################################################################################################################333
  clust = []

  gpus = get_available_gpus()
  gpu_Cnt = len(gpus)
  global ind
  k = max(index1)
  print("Start training and predicting data...")
  for ind in range(0, max(index1), gpu_Cnt):

    with concurrent.futures.ThreadPoolExecutor(len(gpus)) as executor:

      while k > 0:
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
          outputs = pd.concat(outputs, results)
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
  res.to_csv("cnnimpute/output/"+str(args.output))
  stop_end = timeit.default_timer()
  print('total Time: ', stop_end - start_begin)

# if __name__ == "__main__":
#     cnnImpute()