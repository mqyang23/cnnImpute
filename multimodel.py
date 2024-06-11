
import random
import pandas as pd
import numpy as np
from keras import backend as K
from keras.callbacks import CSVLogger
from itertools import chain
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from sklearn.utils import shuffle
from numpy import zeros, newaxis
from keras.models import Sequential,model_from_json
import tensorflow as tf
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers import Dense, Dropout, Flatten, Conv1D, MaxPool1D, Activation, MaxPooling1D, GlobalAveragePooling1D,BatchNormalization
from tensorflow.python.data.util import options as options_lib
from tensorflow.data.experimental import DistributeOptions, AutoShardPolicy




#seperate list of index into chunks
def chunks(l, n):
  # For item i in a range that is a length of l,
  for i in range(0, len(l), n):
    # Create an index range for l of n items:
    yield l[i:i+n]

#RMSE function
def rmse(y_true, y_pred):
      from keras import backend as K
      mask_value=-1
      mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
      y_pred=y_pred*mask
      y_true=y_true*mask
      y_pred=tf.boolean_mask(y_pred,mask)
      y_true=tf.boolean_mask(y_true,mask)
      return K.sqrt(K.mean(K.square(y_pred- y_true), axis=-1))



#mask version for wMSE
def wMSE(y_true, y_pred, binary=False):
    if binary:
        weights = tf.cast(y_true>0, tf.float32)
    else:
        weights = y_true
    return tf.reduce_mean(weights*tf.square(y_true-y_pred))

#MSE with mask
def mse(y_true, y_pred):
      mask_value=-1
      mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
      y_pred=y_pred*mask
      y_true=y_true*mask
      y_pred=tf.boolean_mask(y_pred,mask)
      y_true=tf.boolean_mask(y_true,mask)
      return K.mean(K.square(y_pred- y_true), axis=-1)



#pearson correlation with mask
def pearson_correlation_coefficient(y_true, y_pred):
  epsilon = 10e-5
  mask_value=-1
  mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
  y=y_pred*mask
  x=y_true*mask
  y= tf.boolean_mask(y, mask)
  x= tf.boolean_mask(x, mask)
  mx = K.mean(x)
  my = K.mean(y)
  xm, ym = x - mx, y - my
  r_num = K.sum(xm * ym)
  x_square_sum = K.sum(xm * xm)
  y_square_sum = K.sum(ym * ym)
  r_den = K.sqrt(x_square_sum * y_square_sum)
  r = r_num / (r_den + epsilon)
  return K.mean(r)




#another pearson function
def pearson_r(y_true, y_pred):
  # use smoothing for not resulting in NaN values
  # pearson correlation coefficient
  # https://github.com/WenYanger/Keras_Metrics
  epsilon = 10e-5
  x = y_true
  y = y_pred
  mx = K.mean(x)
  my = K.mean(y)
  xm, ym = x - mx, y - my
  r_num = K.sum(xm * ym)
  x_square_sum = K.sum(xm * xm)
  y_square_sum = K.sum(ym * ym)
  r_den = K.sqrt(x_square_sum * y_square_sum)
  r = r_num / (r_den + epsilon)
  return K.mean(r)




DistributeOptions.auto_shard_policy = options_lib.create_option(
    name="auto_shard_policy",
    ty=AutoShardPolicy,
    docstring="The type of sharding to use. See "
    "`tf.data.experimental.AutoShardPolicy` for additional information.",
    default_factory=lambda: AutoShardPolicy.DATA,
)




def define_model(learning_rate,dropout_rate):

  model = Sequential()

  model.add(Conv1D(filters=16, kernel_size=3, activation='relu', input_shape=(512 * 5, 1)))
  model.add(MaxPooling1D(pool_size=2, strides=2, padding='valid'))
  model.add(BatchNormalization())


  model.add(Conv1D(filters=16, kernel_size=3, activation='relu'))
  model.add(GlobalAveragePooling1D())
  model.add(BatchNormalization())
  model.add(Flatten())
  model.add(Dropout(dropout_rate))
  model.add(Dense(1024, activation="relu", kernel_initializer="glorot_normal"))
  model.add(Dense(512, activation="relu", kernel_initializer="glorot_normal"))

  seed = 1111
  random.seed(seed)
  np.random.seed(seed)
  tf.random.set_seed(seed)
  adam = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, amsgrad=False)

  model.compile(optimizer=adam, loss=mse, metrics=[pearson_correlation_coefficient])

  return model




def loadrowname(ind1):
  ind_list = []
  # read gene name information from the txt file
  with open("/home/junie/cnnImpute/tmp_file/rowname" + str(ind1 + 1) + "_withlabel.txt") as f:
    for tmp_index, line in enumerate(f):
      if ' ' in line:
        ind_list.append(int(line.split('" ')[1]))
  return ind_list

def findinputoutputgenes(ind_list,D,dropoutRate_threshold,inputGene_threshold,geneall0):
  cellimpute = []
  geneselect_database = []

  for nrow in ind_list:  ###axis=1 means row
    if (len(np.where(D[nrow, :] > dropoutRate_threshold)[0]) / D.shape[1] <= inputGene_threshold) and (
            nrow not in geneall0):  ####if the gene have dropout
      geneselect_database.append(nrow)
    if (len(np.where(D[nrow, :] > dropoutRate_threshold)[0]) != 0) and (nrow not in geneall0):
      cellimpute.append(nrow)
  return [cellimpute,geneselect_database]

def inputGene(corr1):
  # find top 5 related gene for each gene
  X = []
  for row in range(0, corr1.shape[0]):
    a1 = corr1.iloc[row, :].sort_values(ascending=False).index.to_numpy()[
         0:5]  ###extract genes that related with
    # a1=[x for x in a1 if x not in cellimpute][0:5] ####extract gene that in a1 but not in cellimpute1 top 5
    X.append(a1)

  # extract all related genes into list
  imputgene_nset = list(chain(*X))
  return imputgene_nset


def runModel(learning_rate,batch_size,max_epochs,dropout_rate,dropoutRate_threshold,inputGene_threshold,data,data1,data_rescaled,ind,geneall0,index1,gpus):
  thre = 0
  # seperate index into different cluster
  for i in range(1, int(max(index1)) + 1):
    text = "cluster" + str(i)  # str means string
    globals()[str(text)] = np.where(index1 == i)[0]

  with tf.compat.v1.Session(graph=tf.Graph()) as sess:

    tf.compat.v1.keras.backend.set_session(sess)

    with tf.device(gpus):

      ind1 = ind + int(gpus[12])
      data_sub = np.array(data[globals()[str("cluster" + str(ind1 + 1))], :])
      # read droprate probability from csv file
      D = pd.read_csv("/home/junie/cnnImpute/tmp_file/droprate"+str(ind1+1)+"_withlabel.csv", index_col=0)
      # convert it into array
      D = np.array(D)

      ind_list = loadrowname(ind1)

      [cellimpute, geneselect_database] = findinputoutputgenes(ind_list,D,dropoutRate_threshold,inputGene_threshold,geneall0)


      # extract this cluster data
      df_tmp = pd.DataFrame(data[globals()[str("cluster" + str(ind1 + 1))], :])
      # extract genes that appears in these two lists from df_tmp
      df_tmp = df_tmp.iloc[:, list(set(geneselect_database) | set(cellimpute))]

      #calculate correlation
      corr = np.corrcoef(np.array(df_tmp).T)
      corr = pd.DataFrame(corr)
      corr.index = df_tmp.columns
      corr.columns = df_tmp.columns

      # seperate datasets into several 50dimensional same size
      cellimpute_nset = list(chunks(cellimpute, 512))

      for nset in range(0, len(cellimpute_nset)):

        if nset == 0:
            globals()["model" + str(int(gpus[12]))] = define_model(learning_rate,dropout_rate)

        cellimpute_sub = cellimpute_nset[nset]

        # find gene that in geneselect_database but not in the cell_impute_sub
        redundance = [int(i) for i in geneselect_database if i not in cellimpute_sub]

        corr1 = corr.loc[cellimpute_sub, redundance]
        corr1 = pd.DataFrame(corr1)
        # extract all related genes into list
        imputgene_nset = inputGene(corr1)
        y_data = data1[globals()[str("cluster" + str(ind1 + 1))], :][:, cellimpute_sub]
        y_data1 = np.copy(y_data)


        subD1 = np.transpose(D)[:, cellimpute_sub]

        for w in range(0, y_data1.shape[1]):
            y_data1[np.where(subD1[:, w] > 0.5)[0], w] = -1



        x_data = data_rescaled.iloc[globals()[str("cluster" + str(ind1 + 1))], imputgene_nset]
        x_data1 = np.copy(x_data)
        # mask data if shape is less than 512
        if (y_data1.shape[1] < 512):
            y_data1 = np.pad(y_data1, ((0, 0), (0, 512 - y_data1.shape[1])), 'constant', constant_values=-1)
            x_data = np.pad(x_data, ((0, 0), (0, 2560 - x_data.shape[1])), 'constant', constant_values=0)
            x_data1 = np.pad(x_data1, ((0, 0), (0, 2560 - x_data1.shape[1])), 'constant', constant_values=0)


        x_data1, y_data1 = shuffle(x_data1, y_data1)
        x_train, x_val, y_train, y_val = train_test_split(x_data1, y_data1, test_size=0.1)
        x_train = np.array(x_train)[:, :, newaxis]
        x_val = np.array(x_val)[:, :, newaxis]
        y_train = np.array(y_train)
        y_val = np.array(y_val)

        seed = 1111
        random.seed(seed)
        np.random.seed(seed)
        tf.random.set_seed(seed)

        csv_logger = CSVLogger("tmp" + str(thre) + '_log.csv', append=True, separator=';')
        monitor_val_acc = EarlyStopping(monitor='val_loss', patience=5, min_delta=thre)

        globals()["model" + str(int(gpus[12]))].fit(x_train, y_train, validation_data=(x_val, y_val), epochs=max_epochs, batch_size=batch_size,
                  callbacks=[csv_logger, monitor_val_acc])



        x_data2 = np.array(x_data)[:, :, newaxis]
        predict = np.array(globals()["model" + str(int(gpus[12]))].predict(x_data2))


        # recover y_data for point we think it is missing(default: dropout probability > 0.5)
        subsetD = np.transpose(D)[:, cellimpute_sub]
        for w in range(0, y_data.shape[1]):
            loca = np.where(subsetD[:, w] > 0.5)[0]

            y_data[loca, w] = predict[loca, w].copy()


        # recover data points in the data
        for ncol in range(0, y_data.shape[1]):
            data_sub[:, cellimpute_sub[ncol]] = y_data[:, ncol].copy()

  return pd.DataFrame(data_sub)

