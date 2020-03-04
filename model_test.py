#!/usr/bin/python3
# -*-coding:utf-8-*-

# @Time     : 2020-03-03 16:21
# @Author   : Sundw
# @File     : model_test_final.py
# @Software : PyCharm
import pandas as pd
import numpy as np
import tensorflow as tf
import scipy.io as scio
import h5py
from tensorflow.contrib import rnn
import argparse
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser()
parser.add_argument("name")
args = parser.parse_args()

# ====================================================================
# 导入数据
k_neighbors = 3
n_class = 2

name = './feature_test/' + args.name
testdata = pd.read_table(name,header=None,sep=" ")
testdata=testdata.drop(labels=0,axis=1)
data_test=testdata.values.astype(np.float32)
nrows = data_test.shape[0]
n = 100 - nrows % 100
zero_data = np.zeros([n,k_neighbors*20+21])
data_test = np.vstack((data_test,zero_data))
m = data_test.shape[0]
test_data = data_test[0:m, 2:20 * k_neighbors + 20]
test_label = data_test[0:m, 20 * k_neighbors + 20]


def convert_to_one_hot(y, n_class):
    y = y.astype(int)
    return np.eye(n_class)[y.reshape(-1)]


test_label = convert_to_one_hot(test_label, n_class)
test_label = test_label.astype(np.float32)
batch_size = 100
label = tf.placeholder(tf.float32, [batch_size, 2])
test_batch = test_label.shape[0] // batch_size

# =======================================================================
# ======== graph convolution =============

filters = 32
k_neighbors = 3
# define center nodes, neighbours, edges
vc_l = tf.placeholder(tf.float32, [batch_size,9])
vc_r = tf.placeholder(tf.float32, [batch_size,9])
vn_l = tf.placeholder(tf.float32, [batch_size,k_neighbors*9])
vn_r = tf.placeholder(tf.float32, [batch_size,k_neighbors*9])
e_l = tf.placeholder(tf.float32, [batch_size, k_neighbors])
e_r = tf.placeholder(tf.float32, [batch_size, k_neighbors])

# test_dict = {label:test_label,vc_l:test_vc_l,vc_r:test_vc_r,vn_l:test_vn_l,vn_r:test_vn_r,e_l:test_e_l,e_r:test_e_r}

# define Wc, Wn, We, b
Wc = tf.Variable(tf.truncated_normal([9,filters], stddev=0.1))
Wn = tf.Variable(tf.truncated_normal([9,filters], stddev=0.1))
We = tf.Variable(tf.truncated_normal([1,filters], stddev=0.1))
b = tf.Variable(tf.truncated_normal([1,filters], stddev=0.1))


def graph_convolution(vc, Wc, vn, Wn, e, We, b, k_neighbors):
    I9 = np.eye(9).astype(np.float32)
    eye1 = tf.concat([I9, I9, I9], 0)
    I3 = tf.constant(1.0, shape=[3, 1], dtype=tf.float32)
    zc = tf.matmul(vc, Wc)
    # z = tf.tanh(z)

    zn = tf.matmul(tf.matmul(vn, eye1), Wn) * 1.0 / k_neighbors

    ze = tf.matmul(tf.matmul(e, I3), We) * 1.0 / k_neighbors

    z = zc + zn + ze + b
    return z


z_l = graph_convolution(vc_l,Wc,vn_l,Wn,e_l,We,b,k_neighbors)
z_r = graph_convolution(vc_r,Wc,vn_r,Wn,e_r,We,b,k_neighbors)
z1 = tf.concat([z_l,z_r],axis = 1)

# ============ mlstm ==============

n_inputs = filters
max_time = 2
lstm_size = 32
filters2 = 16
layer_num = 2
keep_prob=1.0

inputs = tf.reshape(z1,[-1,max_time,n_inputs])
# keep_prob = tf.placeholder(tf.float32)
def unit_lstm():
    # 定义一层LSTM_cell,只需要说明hidden_size,它会自动匹配输入的X的维度
    lstm_cell = rnn.BasicLSTMCell(num_units=lstm_size,forget_bias=1.0,state_is_tuple=True)
    # 添加dropout layer,一般只设置output_keep_prob
    lstm_cell = rnn.DropoutWrapper(cell=lstm_cell,input_keep_prob=keep_prob)
    return lstm_cell
# 调用MultiRNNCell来实现多层LSTM
mlstm_cell = rnn.MultiRNNCell([unit_lstm() for i in range(4)],state_is_tuple=True)

# 用全零来初始化state
init_state = mlstm_cell.zero_state(batch_size, dtype=tf.float32)
outputs, state = tf.nn.dynamic_rnn(mlstm_cell, inputs=inputs,initial_state=init_state, time_major=False)
h_state = outputs[:,-1,:]


weights_lstm = tf.Variable(tf.truncated_normal([lstm_size, filters2], stddev=0.1))
biases_lstm = tf.Variable(tf.constant(0.1, shape=[filters2]))

z2 = tf.nn.tanh(tf.matmul(h_state,weights_lstm) + biases_lstm)

# ========== Dense ============

weights_dense = tf.Variable(tf.truncated_normal([filters2,n_class], stddev=0.1))
biases_dense = tf.Variable(tf.constant(1.0, shape=[2]))

prediction = tf.nn.softmax(tf.matmul(z2,weights_dense) + biases_dense)

# ========== Loss and optimizer ===========

# loss function
cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=label, logits = prediction))

# AdamOptimizer优化
train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)

correct_prediction = tf.equal(tf.argmax(label,1), tf.argmax(prediction,1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

saver = tf.train.Saver()

with tf.Session() as sess:
    # sess.run(init)

    saver.restore(sess, "./model_341/model.ckpt")

    test_data_all = test_data.copy()
    test_label_all = test_label.copy()
    pred = np.array([[0, 0]], dtype=np.float32)

    for batch in range(test_batch):
        test_data_batch = test_data_all[0:batch_size, :]
        test_label_batch = test_label_all[0:batch_size, :]
        if batch != test_batch - 1:
            test_data_all = test_data_all[batch_size:, :]
            test_label_all = test_label_all[batch_size:, :]

        test_vc_l = test_data_batch[:, 0:9]
        test_vn_l = test_data_batch[:, 9:9 * k_neighbors + 9]
        test_e_l = test_data_batch[:, 9 * k_neighbors + 9:10 * k_neighbors + 9]
        test_vc_r = test_data_batch[:, 10 * k_neighbors + 9:10 * k_neighbors + 18]
        test_vn_r = test_data_batch[:, 10 * k_neighbors + 18:19 * k_neighbors + 18]
        test_e_r = test_data_batch[:, 19 * k_neighbors + 18:20 * k_neighbors + 18]

        test_dict = {label: test_label_batch, vc_l: test_vc_l, vc_r: test_vc_r, vn_l: test_vn_l, vn_r: test_vn_r,
                     e_l: test_e_l, e_r: test_e_r}
        temp = sess.run(prediction, feed_dict=test_dict)
        pred = np.concatenate((pred, temp), axis=0)

    pred1 = pred[1:nrows+1]
    # pred = sess.run(prediction, feed_dict=test_dict
    # np.savetxt('./results/1bv4_AB.txt', pred1, fmt="%s", delimiter=" ")


print('The program has been finished.')

resnum = data_test[0:nrows,0:2]
data_label = data_test[0:nrows,-1]
datalabel = data_label.reshape(-1,1)
result = np.concatenate((pred1,resnum,datalabel),axis=1)
result = result[np.argsort(result[:,0])]
result = result[:,2:]

savename = './results/' + args.name
np.savetxt(savename, result, fmt="%s", delimiter=" ")