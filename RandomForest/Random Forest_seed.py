import pandas as pd
import numpy as np
np.random.seed(1980) # S'ajusta llavor per tal que el resultat sigui reproduible
from keras.models import load_model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import f1_score

# Es carrega el model d'encoder. El model final no 
# utilitza l'autoencoder, es mante per a futures iteracions
encoder = load_model('autoencoder/encoder_H20.h5')

# Es carreguen les dades train i es passen per l'encoder. 
train = pd.read_csv('data/trainPredictNorm.csv')
train = np.array(train)

train_onehot = train[:, 11:]
##encoded_train = encoder.predict(train_onehot)

# S'afegeix la informacio referent a la llavor
encoded_train = np.column_stack((train[:, 6:11], train_onehot))

# Es carreguen les etiquetes del subset train
y_train = train[:, 0]
y_train = y_train.astype("float64")

# Es carreguen les dades test i es passen per l'encoder
test = pd.read_csv('data/testPredictNorm.csv')
test = np.array(test)

test_onehot = test[:, 11:]
##encoded_test = encoder.predict(test_onehot)

# S'afegeix la informacio referent a la llavor
encoded_test = np.column_stack((test[:, 6:11], test_onehot))

# Es crea el model amb Random Forest
model= RandomForestClassifier(n_estimators=35, bootstrap=False)

# S'entrena el model amb el subset train
model.fit(encoded_train, y_train)

# Es guarda el model per ser utilitzar posteriorment
import pickle
filename = 'models/RF_model.sav'
pickle.dump(model, open(filename, 'wb'))

# Estimem la precisio del model per cross-validation
scores = cross_val_score(model, encoded_train, y_train, cv=10)
print scores
print("Accuracy estimation: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))

# Es prediuen les categories dels exemples del subset test i es guarden per ser avaluades
predicted = model.predict(encoded_test)

np.savetxt("results/RF_seed.csv", predicted, delimiter=",")

# Es calcula l'exactitud i el F1-score del model
y_test = test[:, 0]
y_test = y_test.astype("float64")
result = model.score(encoded_test, y_test)
print("Accuracy: %0.4f" % (result))

F1 = f1_score(y_test, predicted, pos_label=1)
print F1
