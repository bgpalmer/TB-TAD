import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, mean_absolute_error, mean_squared_error 
from sklearn.model_selection import train_test_split

from sklearn.ensemble import GradientBoostingClassifier
from model import PytorchBasedGenericGradientBoost


# Definition of Hyper-Parameters
NUM_CLASSIFIERS = 2
MAX_DEPTH = 10
GRADIENT_BOOST_LEARNING_RATE = 0.1
MINIMIZER_LEARNING_RATE = 0.005
MINIMIZER_TRAINING_EPOCHS = 10000

# Read the training data


# Change this again for Live Data
df = pd.read_csv("./all_data.csv", sep=",")

# combination on me
me_cols = ['me','me_L_1','me_L_2','me_L_3','me_L_4','me_L_5','me_L_6','me_L_7','me_L_8','me_L_9','me_L_10','me_R_1','me_R_2','me_R_3','me_R_4','me_R_5','me_R_6','me_R_7','me_R_8','me_R_9','me_R_10']
chr = df.loc[:, ['chr']]
new_df = pd.DataFrame()
new_df['me'] = df[me_cols].values.tolist()
new_df['me'] = new_df['me'].apply(lambda x: np.array(x, dtype=float))

# Bad Solution 
new_df['me_bad_solution'] = new_df['me'].apply(lambda x: np.mean(x))
###

new_df['chr'] = chr

# X = new_df.loc[:,['me','chr']]
X = new_df.loc[:,['me_bad_solution']]
y = df.loc[:, "bds"]
X = np.nan_to_num(X, 0)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

print(np.array(X_train))

#print(np.asarray(np.array(X_train)))


# Running the custom algorithm 
# custom = PytorchBasedGenericGradientBoost("classifier", NUM_CLASSIFIERS, MAX_DEPTH, GRADIENT_BOOST_LEARNING_RATE=GRADIENT_BOOST_LEARNING_RATE, MINIMIZER_LEARNING_RATE=MINIMIZER_LEARNING_RATE, MINIMIZER_TRAINING_EPOCHS=MINIMIZER_TRAINING_EPOCHS)
# custom.fit(X_train, y_train)
# predictions_train = custom.predict(X_train)
# predictions_test = custom.predict(X_test)
# print("Custom Implementation : Accuracy score for training data : {}".format(accuracy_score(np.round(predictions_train), y_train)))
# print("Custom Implementation : Accuracy score for testing data : {}".format(accuracy_score(np.round(predictions_test), y_test)))

# Running the vanilla sklearn algorithm
classifier = GradientBoostingClassifier(n_estimators=NUM_CLASSIFIERS, max_depth=MAX_DEPTH)
classifier.fit(np.array(X_train), np.array(y_train))
y_pred_gb = classifier.predict(X_test)
print("Vanilla Implementation : Accuracy score for training data : {}".format(accuracy_score(classifier.predict(X_train), y_train)))
print("Vanilla Implementation : Accuracy score for testing data : {}".format(accuracy_score(classifier.predict(X_test), y_test)))

mae_gb = mean_absolute_error(y_test, y_pred_gb)

mse_gb = mean_squared_error(y_test, y_pred_gb)

print(f'mae_gb: {mae_gb}, \nmse_gb: {mse_gb}')