from __future__ import division
from __future__ import print_function
from operator import itemgetter
from itertools import combinations
import time
import os
import logging
import sys

sys.path.insert(0, './')

import tensorflow as tf
import numpy as np
import networkx as nx
import scipy.sparse as sp
from sklearn import metrics

from decagon.deep.optimizer import DecagonOptimizer
from decagon.deep.model import DecagonModel
from decagon.deep.minibatch import EdgeMinibatchIterator
from decagon.utility import rank_metrics, preprocessing
from polypharmacy.utility import *
from collections import Counter

from decagon.utility.visualization import WriterTensorboardX
import datetime

# Train on CPU (hide GPU) due to memory constraints
# os.environ['CUDA_VISIBLE_DEVICES'] = ""

# Train on GPU
os.environ["CUDA_DEVICE_ORDER"] = 'PCI_BUS_ID'
os.environ["CUDA_VISIBLE_DEVICES"] = '0'
config = tf.ConfigProto()
config.gpu_options.allow_growth = True

np.random.seed(0)


###########################################################
#
# Functions
#
###########################################################


def get_accuracy_scores(edges_pos, edges_neg, edge_type):
    feed_dict.update({placeholders['dropout']: 0})
    feed_dict.update({placeholders['batch_edge_type_idx']: minibatch.edge_type2idx[edge_type]})
    feed_dict.update({placeholders['batch_row_edge_type']: edge_type[0]})
    feed_dict.update({placeholders['batch_col_edge_type']: edge_type[1]})
    rec = sess.run(opt.predictions, feed_dict=feed_dict)

    def sigmoid(x):
        return 1. / (1 + np.exp(-x))

    # Predict on test set of edges
    preds = []
    actual = []
    predicted = []
    edge_ind = 0
    for u, v in edges_pos[edge_type[:2]][edge_type[2]]:
        score = sigmoid(rec[u, v])
        preds.append(score)
        assert adj_mats_orig[edge_type[:2]][edge_type[2]][u, v] == 1, 'Problem 1'

        actual.append(edge_ind)
        predicted.append((score, edge_ind))
        edge_ind += 1

    preds_neg = []
    for u, v in edges_neg[edge_type[:2]][edge_type[2]]:
        score = sigmoid(rec[u, v])
        preds_neg.append(score)
        assert adj_mats_orig[edge_type[:2]][edge_type[2]][u, v] == 0, 'Problem 0'

        predicted.append((score, edge_ind))
        edge_ind += 1

    preds_all = np.hstack([preds, preds_neg])
    preds_all = np.nan_to_num(preds_all)
    labels_all = np.hstack([np.ones(len(preds)), np.zeros(len(preds_neg))])
    predicted = list(zip(*sorted(predicted, reverse=True, key=itemgetter(0))))[1]

    roc_sc = metrics.roc_auc_score(labels_all, preds_all)
    aupr_sc = metrics.average_precision_score(labels_all, preds_all)
    apk_sc = rank_metrics.apk(actual, predicted, k=50)

    return roc_sc, aupr_sc, apk_sc


def construct_placeholders(edge_types):
    placeholders = {
        'batch': tf.placeholder(tf.int32, name='batch'),
        'batch_edge_type_idx': tf.placeholder(tf.int32, shape=(), name='batch_edge_type_idx'),
        'batch_row_edge_type': tf.placeholder(tf.int32, shape=(), name='batch_row_edge_type'),
        'batch_col_edge_type': tf.placeholder(tf.int32, shape=(), name='batch_col_edge_type'),
        'degrees': tf.placeholder(tf.int32),
        'dropout': tf.placeholder_with_default(0., shape=()),
    }
    placeholders.update({
        'adj_mats_%d,%d,%d' % (i, j, k): tf.sparse_placeholder(tf.float32)
        for i, j in edge_types for k in range(edge_types[i, j])})
    placeholders.update({
        'feat_%d' % i: tf.sparse_placeholder(tf.float32)
        for i, _ in edge_types})
    return placeholders


###########################################################
#
# Load and preprocess data (This is a dummy toy example!)
#
###########################################################

####
# The following code uses artificially generated and very small networks.
# Expect less than excellent performance as these random networks do not have any interesting structure.
# The purpose of main.py is to show how to use the code!
#
# All preprocessed datasets used in the drug combination study are at: http://snap.stanford.edu/decagon:
# (1) Download datasets from http://snap.stanford.edu/decagon to your local machine.
# (2) Replace dummy toy datasets used here with the actual datasets you just downloaded.
# (3) Train & test the model.
####

combo2stitch, combo2se, se2name = load_combo_se(fname='./data/csv/bio-decagon-combo.csv')
gene_net, gene_idx_dict = load_ppi('./data/csv/bio-decagon-ppi.csv')
stitch2proteins = load_targets(fname='./data/csv/bio-decagon-targets.csv')
stitch2se, se2name_mono = load_mono_se('./data/csv/bio-decagon-mono.csv')

all_drug = []
for drug_pair in combo2stitch.values():
    all_drug += drug_pair

# look-up table for drugs
drug_id = set(all_drug)
print('Drug number = %d' % len(drug_id))
drug_idx_dict = {item: idx for idx, item in enumerate(drug_id)}


# we only focus on most common 964 side effect
def get_se_counter(se_map):
    side_effects = []
    for drug in se_map:
        side_effects += list(set(se_map[drug]))
    return Counter(side_effects)


combo_counter = get_se_counter(combo2se)

common_se_combo_id = []
common_se_combo_id_counts = []
common_se_combo_names = []
for se, count in combo_counter.most_common(964):
    common_se_combo_id += [se]
    common_se_combo_id_counts += [count]
    common_se_combo_names += [se2name[se]]

# look-up table for target combo side effect
se_id = common_se_combo_id
se_idx_dict = {item: idx for idx, item in enumerate(se_id)}

val_test_size = 0.1
n_genes = len(gene_idx_dict)
n_drugs = len(drug_id)
n_drugdrug_rel_types = len(se_id)

gene_adj = nx.adjacency_matrix(gene_net)
gene_degrees = np.array(gene_adj.sum(axis=0)).squeeze()

gene_row = []
drug_col = []

for drug_selected, target_gene in stitch2proteins.items():
    drug_selected_idx = drug_idx_dict[drug_selected]

    for gene in target_gene:
        try:
            target_gene_idx = gene_idx_dict[gene]
            gene_row += [target_gene_idx]
            drug_col += [drug_selected_idx]
        except:
            # some target proteins lie outside ppi graph
            pass

row = np.array(gene_row)
col = np.array(drug_col)
data = np.ones_like(row)
gene_drug_adj = sp.csr_matrix((data, (row, col)), shape=(len(gene_idx_dict), len(drug_id)))
drug_gene_adj = gene_drug_adj.transpose(copy=True)

drug_drug_adj_list = np.zeros([n_drugdrug_rel_types, n_drugs, n_drugs])
for drug_pair, se_list in combo2se.items():
    drug_1, drug_2 = combo2stitch[drug_pair]
    drug_1_id, drug_2_id = drug_idx_dict[drug_1], drug_idx_dict[drug_2]

    for se in se_list:
        if se in se_idx_dict:
            se_idx = se_idx_dict[se]
            drug_drug_adj_list[se_idx][drug_1_id, drug_2_id] = 1
            drug_drug_adj_list[se_idx][drug_2_id, drug_1_id] = 1

drug_drug_adj_list = [sp.csr_matrix(mat) for mat in drug_drug_adj_list]
drug_degrees_list = [np.array(drug_adj.sum(axis=0)).squeeze() for drug_adj in drug_drug_adj_list]

# data representation
adj_mats_orig = {
    (0, 0): [gene_adj],
    (0, 1): [gene_drug_adj],
    (1, 0): [drug_gene_adj],
    (1, 1): drug_drug_adj_list,
}
degrees = {
    0: [gene_degrees],
    1: drug_degrees_list,
}

# featureless (genes)
gene_feat = sp.identity(n_genes)
gene_nonzero_feat, gene_num_feat = gene_feat.shape
gene_feat = preprocessing.sparse_to_tuple(gene_feat.tocoo())

# features (drugs): side effects of individual drugs were used as additional features for drug nodes.
mono_se_list = []
for mono_se in stitch2se.values():
    mono_se_list.append(list(mono_se))

mono_se_list = np.concatenate(mono_se_list)

mono_se_id = np.array(list(set(mono_se_list)))
mono_se_idx_dict = {item: idx for idx, item in enumerate(mono_se_id)}

drug_feature_row = []
drug_feature_col = []

for drug_selected, mono_se_list in stitch2se.items():
    drug_selected_idx = drug_idx_dict[drug_selected]

    for mono_se in mono_se_list:
        mono_se_idx = mono_se_idx_dict[mono_se]
        drug_feature_row += [drug_selected_idx]
        drug_feature_col += [mono_se_idx]

row = np.array(drug_feature_row)
col = np.array(drug_feature_col)
data = np.ones_like(row)
drug_feat = sp.csr_matrix((data, (row, col)), shape=(n_drugs, len(mono_se_id)))

drug_nonzero_feat, drug_num_feat = drug_feat.shape
drug_feat = preprocessing.sparse_to_tuple(drug_feat.tocoo())

# data representation
num_feat = {  # input dim
    0: gene_num_feat,
    1: drug_num_feat,
}
nonzero_feat = {
    0: gene_nonzero_feat,
    1: drug_nonzero_feat,
}
feat = {
    0: gene_feat,
    1: drug_feat,
}

edge_type2dim = {k: [adj.shape for adj in adjs] for k, adjs in adj_mats_orig.items()}
edge_type2decoder = {
    (0, 0): 'bilinear',
    (0, 1): 'bilinear',
    (1, 0): 'bilinear',
    (1, 1): 'dedicom',  # To capture the polypharmacy combinatorics
}

edge_types = {k: len(v) for k, v in adj_mats_orig.items()}
num_edge_types = sum(edge_types.values())
print("Edge types:", "%d" % num_edge_types)

###########################################################
#
# Settings and placeholders
#
###########################################################

flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_integer('neg_sample_size', 1, 'Negative sample size.')
flags.DEFINE_float('learning_rate', 0.001, 'Initial learning rate.')
flags.DEFINE_integer('epochs', 100, 'Number of epochs to train.')
flags.DEFINE_integer('hidden1', 64, 'Number of units in hidden layer 1.')
flags.DEFINE_integer('hidden2', 32, 'Number of units in hidden layer 2.')
flags.DEFINE_float('weight_decay', 0, 'Weight for L2 loss on embedding matrix.')
flags.DEFINE_float('dropout', 0.1, 'Dropout rate (1 - keep probability).')
flags.DEFINE_float('max_margin', 0.1, 'Max margin parameter in hinge loss')
flags.DEFINE_integer('batch_size', 512, 'minibatch size.')
flags.DEFINE_boolean('bias', True, 'Bias term.')
# Important -- Do not evaluate/print validation performance every iteration as it can take
# substantial amount of time
PRINT_PROGRESS_EVERY = 500

print("Defining placeholders")
placeholders = construct_placeholders(edge_types)

###########################################################
#
# Create minibatch iterator, model and optimizer
#
###########################################################

print("Create minibatch iterator")
minibatch = EdgeMinibatchIterator(
    adj_mats=adj_mats_orig,
    feat=feat,
    edge_types=edge_types,
    batch_size=FLAGS.batch_size,
    val_test_size=val_test_size
)

print("Create model")
model = DecagonModel(
    placeholders=placeholders,
    num_feat=num_feat,
    nonzero_feat=nonzero_feat,
    edge_types=edge_types,
    decoders=edge_type2decoder,
)

print("Create optimizer")
with tf.name_scope('optimizer'):
    opt = DecagonOptimizer(
        embeddings=model.embeddings,
        latent_inters=model.latent_inters,
        latent_varies=model.latent_varies,
        degrees=degrees,
        edge_types=edge_types,
        edge_type2dim=edge_type2dim,
        placeholders=placeholders,
        batch_size=FLAGS.batch_size,
        margin=FLAGS.max_margin
    )

print("Initialize session")
sess = tf.Session()
sess.run(tf.global_variables_initializer())
feed_dict = {}

# setup visualization writer instance

start_time = datetime.datetime.now().strftime('%m%d_%H%M%S')
writer_dir = os.path.join('./experiment/log/', start_time)

if not os.path.exists(writer_dir):
    os.makedirs(writer_dir)

writer = WriterTensorboardX(writer_dir, logging.getLogger(start_time), enable=True)


###########################################################
#
# epoch evaluation
#
###########################################################

def epoch_eval(minibatch, epoch, mode='val'):
    # test #
    test_list = [1, 2, 3, 4, 5]
    writer.add_scalars('ave_test', np.mean(test_list))
    writer.add_histogram('distribution_test', np.array(test_list))

    side_effect_roc_score = []
    side_effect_auprc_score = []
    side_effect_apk_score = []
    for et in range(num_edge_types):
        if mode == 'val':
            roc_score, auprc_score, apk_score = get_accuracy_scores(
                minibatch.val_edges, minibatch.val_edges_false, minibatch.idx2edge_type[et])
        if mode == 'test':
            roc_score, auprc_score, apk_score = get_accuracy_scores(
                minibatch.test_edges, minibatch.test_edges_false, minibatch.idx2edge_type[et])
        #     print("Edge type=", "[%02d, %02d, %02d]" % minibatch.idx2edge_type[et])
        #     print("Edge type:", "%04d" % et, "Test AUROC score", "{:.5f}".format(roc_score))
        #     print("Edge type:", "%04d" % et, "Test AUPRC score", "{:.5f}".format(auprc_score))
        #     print("Edge type:", "%04d" % et, "Test AP@k score", "{:.5f}".format(apk_score))

        writer.set_step(epoch, mode)
        if et <= 2:
            writer.add_scalars('%04d_roc_score' % (et), roc_score)
            writer.add_scalars('%04d_auprc_score' % (et), auprc_score)
            writer.add_scalars('%04d_apk_score' % (et), apk_score)
        else:
            side_effect_roc_score.append(roc_score)
            side_effect_auprc_score.append(auprc_score)
            side_effect_apk_score.append(apk_score)

    writer.add_scalars('ave_side_effect_auc_score', np.mean(side_effect_roc_score))
    writer.add_scalars('ave_side_effect_auprc_score', np.mean(side_effect_auprc_score))
    writer.add_scalars('ave_side_effect_apk_score', np.mean(side_effect_apk_score))

    writer.add_histogram('distribution_side_effect_auc_score', np.array(side_effect_roc_score))
    writer.add_histogram('distribution_side_effect_auprc_score', np.array(side_effect_auprc_score))
    writer.add_histogram('distribution_side_effect_apk_score', np.array(side_effect_apk_score))


###########################################################
#
# Train model
#
###########################################################

print("Train model")
for epoch in range(FLAGS.epochs):

    minibatch.shuffle()

    itr = 0
    while not minibatch.end():
        # Construct feed dictionary
        feed_dict = minibatch.next_minibatch_feed_dict(placeholders=placeholders)
        feed_dict = minibatch.update_feed_dict(
            feed_dict=feed_dict,
            dropout=FLAGS.dropout,
            placeholders=placeholders)

        t = time.time()

        # Training step: run single weight update
        outs = sess.run([opt.opt_op, opt.cost, opt.batch_edge_type_idx], feed_dict=feed_dict)
        train_cost = outs[1]
        batch_edge_type = outs[2]

        if itr % PRINT_PROGRESS_EVERY == 0:
            val_auc, val_auprc, val_apk = get_accuracy_scores(
                minibatch.val_edges, minibatch.val_edges_false,
                minibatch.idx2edge_type[minibatch.current_edge_type_idx])

            print("Epoch:", "%04d" % (epoch + 1), "Iter:", "%04d" % (itr + 1), "Edge:", "%04d" % batch_edge_type,
                  "train_loss=", "{:.5f}".format(train_cost),
                  "val_roc=", "{:.5f}".format(val_auc), "val_auprc=", "{:.5f}".format(val_auprc),
                  "val_apk=", "{:.5f}".format(val_apk), "time=", "{:.5f}".format(time.time() - t))

        if itr == 0:  # test model for each epoch
            print("test model !!!")
            epoch_eval(minibatch, epoch, mode='test')

        itr += 1

print("Optimization finished!")
