import deepchem as dc
from deepchem import splits
import tensorflow.keras as keras
import matplotlib.pyplot as plt
import logging

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")

logging.info("Opening dataset and splitting")
dataset = dc.data.DiskDataset("data/dataset")
splitter = splits.RandomSplitter()
train, test = splitter.train_test_split(dataset, seed=42)

logging.info("Initializing Model")
model = tf.keras.Sequential([
    keras.layers.Dense(units=32, activation='relu', input_shape=dc_dataset.X.shape[1:]),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(units=32, activation='relu'), 
    keras.layers.Dropout(0.2),
    keras.layers.Dense(units=1),
])

model.compile(loss='mae', optimizer='adam')
print(model.summary())

history = model.fit(
    train.X, train.y,
    validation_data=(test.X,test.y),
    batch_size=100,
    epochs=30
)

history_df = pd.DataFrame(history.history)
history_df[['loss', 'val_loss']].plot()
