import keras
from keras import Input
from keras.engine import Model
from keras.layers import Bidirectional, GRU, Dense, Dropout, LSTM, SimpleRNN
from keras import metrics
from keras.layers.advanced_activations import PReLU, LeakyReLU

import train.CustomLoss as cl

METRICS = ['mae', metrics.binary_accuracy]

def simple_model(shape1, shape2, output_activation='sigmoid'):
    """
    Simple Model architecture
    """
    # PEPTIDE
    pept_input = Input(shape=shape1, name='pept_input')
    pept_layer = LSTM(512, dropout=0.1, recurrent_dropout=0.1, name='pept_lstm_1')(pept_input)
    pept_layer = LeakyReLU()(pept_layer)

    # MHC
    mhc_input = Input(shape=shape2, name='mhc_input')
    mhc_layer = LSTM(512, dropout=0.1, recurrent_dropout=0.1, name='mhc_lstm_1')(mhc_input)
    mhc_layer = LeakyReLU()(mhc_layer)

    # COMBINED
    merged = keras.layers.concatenate([pept_layer, mhc_layer])
    merged = Dense(512)(merged)
    merged = LeakyReLU()(merged)
    merged = Dropout(0.4)(merged)
    main_output = Dense(1, activation=output_activation, name='main_output')(merged)

    # GENERATE MODEL
    model = Model(inputs=[mhc_input, pept_input], outputs=[main_output])

    model.compile(optimizer='adam',
                  loss=cl.weighted_binary_loss,
                  metrics=METRICS)
    return model


def model_narrow_deep(shape1, shape2, output_activation='sigmoid'):
    """
    Deep Narrow architecture
    """
    # PEPTIDE
    pept_input = Input(shape=shape1, name='pept_input')

    pept_layer = LSTM(256, name='pept_lstm_1', dropout=0.1, recurrent_dropout=0.1)(pept_input)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dense(256)(pept_layer)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dropout(0.4)(pept_layer)
    pept_layer = Dense(128)(pept_layer)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dropout(0.4)(pept_layer)
    pept_layer = Dense(64)(pept_layer)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dropout(0.4)(pept_layer)

    # MHC
    mhc_input = Input(shape=shape2, name='mhc_input')

    mhc_layer = LSTM(256, name='mhc_lstm_1', dropout=0.1, recurrent_dropout=0.1)(mhc_input)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dense(256)(mhc_layer)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dropout(0.4)(mhc_layer)
    mhc_layer = Dense(128)(mhc_layer)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dropout(0.4)(mhc_layer)
    mhc_layer = Dense(64)(mhc_layer)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dropout(0.4)(mhc_layer)

    # COMBINED
    merged = keras.layers.concatenate([pept_layer, mhc_layer])

    merged = Dense(128)(merged)
    merged = LeakyReLU()(merged)
    merged = Dropout(0.4)(merged)
    merged = Dense(64)(merged)
    merged = LeakyReLU()(merged)
    merged = Dropout(0.4)(merged)
    merged = Dense(32)(merged)
    merged = LeakyReLU()(merged)
    merged = Dropout(0.4)(merged)

    # Finally we add the main logistic regression layer
    main_output = Dense(1, activation=output_activation, name='main_output')(merged)

    # Generate model
    model = Model(inputs=[mhc_input, pept_input], outputs=[main_output])

    model.compile(optimizer='adam',
                  loss=cl.weighted_binary_loss,
                  metrics=METRICS)

    return model


def model_wide_shallow(shape1, shape2, output_activation='sigmoid'):
    """
    Wide Shallow architecture
    """
    # PEPTIDE
    pept_input = Input(shape=shape1, name='pept_input')

    pept_layer = LSTM(512, name='pept_lstm_1', dropout=0.1, recurrent_dropout=0.1)(pept_input)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dense(256)(pept_layer)
    pept_layer = LeakyReLU()(pept_layer)
    pept_layer = Dropout(0.4)(pept_layer)

    # MHC
    mhc_input = Input(shape=shape2, name='mhc_input')

    mhc_layer = LSTM(512, name='mhc_lstm_1', dropout=0.1, recurrent_dropout=0.1)(mhc_input)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dense(256)(mhc_layer)
    mhc_layer = LeakyReLU()(mhc_layer)
    mhc_layer = Dropout(0.4)(mhc_layer)

    merged = keras.layers.concatenate([pept_layer, mhc_layer])
    # We stack a deep densely-connected network on top
    merged = Dense(256)(merged)
    merged = LeakyReLU()(merged)
    merged = Dropout(0.4)(merged)

    # And finally we add the main logistic regression layer
    main_output = Dense(1, activation=output_activation, name='main_output')(merged)

    # Generate Model
    model = Model(inputs=[mhc_input, pept_input], outputs=[main_output])
    model.compile(optimizer='adam',
                  loss=cl.weighted_binary_loss,
                  metrics=METRICS)

    return model
