import tensorflow as tf
import networkx as nx
import numpy as np

def GW_cut(G, root=None, tol=1.E-6):
    @tf.function
    def GW_objective(adj, w, x):
        x = tf.math.l2_normalize(x, axis=1)
        return tf.math.reduce_sum((adj*tf.math.sigmoid(w)) * (1-x @ tf.keras.backend.transpose(x)))
    adj = nx.to_numpy_array(G)
    k = len(G)
    x = tf.Variable(np.random.random((len(G), k)))
    w = tf.Variable(np.random.random((len(G), len(G))))
    optimizer = tf.optimizers.SGD(learning_rate=0.1)
    loss_prev = 0
    for iters in range(100000):
        with tf.GradientTape() as tape:
            res = -GW_objective(adj, w, x)
        gradients = tape.gradient(res, [w, x])
        optimizer.apply_gradients(zip(gradients, [w, x]))
        loss_value = -res.numpy()
        if loss_value-loss_prev < tol:
            break
        loss_prev = loss_value
    u = tf.math.l2_normalize(x, axis=1)
    if root is None:
        r = np.random.random(k)
        r /= np.dot(r,r)
    else:
        r = u[list(G).index(root)]
    return [node for i, node in enumerate(list(G)) if np.dot(r,u[i])>=0]