import tensorflow as tf
import keras.backend as K

import math

from tensorflow.python.framework import constant_op
from tensorflow.python.framework import dtypes
from tensorflow.python.framework import ops
from tensorflow.python.ops import array_ops
from tensorflow.python.ops import candidate_sampling_ops
from tensorflow.python.ops import embedding_ops
from tensorflow.python.ops import gen_nn_ops
from tensorflow.python.ops import math_ops
from tensorflow.python.ops import nn_ops
from tensorflow.python.ops import sparse_ops
from tensorflow.python.ops import variables


_EPSILON = 10e-8


def _to_tensor(x, dtype):
    """Convert the input `x` to a tensor of type `dtype`.

    # Arguments
        x: An object to be converted (numpy array, list, tensors).
        dtype: The destination type.

    # Returns
        A tensor.
    """
    x = tf.convert_to_tensor(x)
    if x.dtype != dtype:
        x = tf.cast(x, dtype)
    return x


def __binary_crossentropy(output, target, pos_weight, from_logits=False):
    """Binary crossentropy between an output tensor and a target tensor.

    # Arguments
        output: A tensor.
        target: A tensor with the same shape as `output`.
        from_logits: Whether `output` is expected to be a logits tensor.
            By default, we consider that `output`
            encodes a probability distribution.

    # Returns
        A tensor.
    """
    # Note: tf.nn.sigmoid_cross_entropy_with_logits
    # expects logits, Keras expects probabilities.
    if not from_logits:
        # transform back to logits
        epsilon = _to_tensor(_EPSILON, output.dtype.base_dtype)
        output = tf.clip_by_value(output, epsilon, 1 - epsilon)
        output = tf.log(output / (1 - output))

    return tf.nn.weighted_cross_entropy_with_logits(targets=target,
                                                    logits=output,
                                                    pos_weight=pos_weight)


def weighted_binary_loss(y_true, y_pred, pos_weight=1/2):
    """
    Custom weighted binary loss function. Penalize negatives 2 times harder than positives.
    """
    # Pos_weight has to be < 1 to penalize false positives!
    return K.mean(__binary_crossentropy(y_pred, y_true, pos_weight=pos_weight), axis=-1)

