import pandas as pd
import numpy as np
np.random.seed(1980) # S'ajusta llavor per tal que el resultat sigui reproduible
from keras.layers import Input, Dense
from keras.models import Model
from keras import regularizers
from keras.callbacks import EarlyStopping
import matplotlib.pyplot as plt

# S'estableix que l'aprenentatge s'aturi si en 10 epoques el model no millora
early_stopping = EarlyStopping(monitor='val_loss', patience=10)

# Es carreguen les dades per entrenar i testar el model
train = pd.read_csv('Data/trainFeatures.csv')
train = np.array(train)
test = pd.read_csv('Data/testFeatures.csv')
test = np.array(test)

# Es construeix l'stacked autoencoder
# funcions d'activacio: 'relu', 'softmax', 'elu', 'softplus', 'softsign', 'tanh', 'sigmoid',
# 'hard_sigmoid', 'linear'
inputs = Input(shape=(220,))
h = Dense(130, activation='linear')(inputs)
encoded = Dense(60, activation='linear', kernel_regularizer=regularizers.l2(1e-5))(h)
h = Dense(130, activation='sigmoid')(encoded)
outputs = Dense(220, activation='sigmoid')(h)

# Es defineix el model, aixi com la funcio 'loss' que permetra l'aprenentatge
encoder = Model(inputs, encoded)
autoencoder = Model(input=inputs, output=outputs)
autoencoder.compile(optimizer='Nadam', loss='mean_squared_error')

# S'entrena el model
autoencoder.fit(train, train,
                epochs=1000,
                shuffle=True,
                validation_data=(test, test),
                callbacks=[early_stopping])

encoder.save('Models/S_encoder_H130_H60.h5')

# serialize model to YAML
S_encoder_H130_H60 = encoder.to_yaml()
with open("Models/S_encoder_H130_H60.yaml", "w") as yaml_file:
    yaml_file.write(S_encoder_H130_H60)

# Es fan les prediccions i es guarden en un arxiu .csv per ser evaluades posteriorment
predictions = autoencoder.predict(test)
np.savetxt("Results/Hidden130-60/decoded_linear.csv", predictions, delimiter=",")

# Es redueixen dimensionalment les dades
encoded_data = encoder.predict(test)
print encoded_data.shape

# Es representa graficament l'evolucio de l'error durant l'entrenament
plt.figure(figsize=(6, 3))
plt.plot(autoencoder.history.history['loss'])
plt.ylabel('error')
plt.xlabel('iteration')
plt.title('training error')
plt.show()