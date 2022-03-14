import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, mean_absolute_error, mean_squared_error 
from sklearn.model_selection import train_test_split
import torch
from sklearn.ensemble import GradientBoostingClassifier
from model import PytorchBasedGenericGradientBoost


# Definition of Hyper-Parameters
NUM_CLASSIFIERS = 2
MAX_DEPTH = 10
GRADIENT_BOOST_LEARNING_RATE = 0.1
MINIMIZER_LEARNING_RATE = 0.005
MINIMIZER_TRAINING_EPOCHS = 10000

# Load in data
df = pd.read_csv("./all_data.csv", sep=",")
new_df = pd.DataFrame()
# Combination Lists
me_cols = ['me','me_L_1','me_L_2','me_L_3','me_L_4','me_L_5','me_L_6','me_L_7','me_L_8','me_L_9','me_L_10','me_R_1','me_R_2','me_R_3','me_R_4','me_R_5','me_R_6','me_R_7','me_R_8','me_R_9','me_R_10']
hist_cols = ['H3K27me3_10A','H3K27me3_10A_L_1','H3K27me3_10A_L_2','H3K27me3_10A_L_3','H3K27me3_10A_L_4','H3K27me3_10A_L_5','H3K27me3_10A_L_6','H3K27me3_10A_L_7','H3K27me3_10A_L_8','H3K27me3_10A_L_9','H3K27me3_10A_L_10','H3K27me3_10A_R_1','H3K27me3_10A_R_2','H3K27me3_10A_R_3','H3K27me3_10A_R_4','H3K27me3_10A_R_5','H3K27me3_10A_R_6','H3K27me3_10A_R_7','H3K27me3_10A_R_8','H3K27me3_10A_R_9','H3K27me3_10A_R_10']
ctcf_cols = ['CTCF_10A','CTCF_10A_L_1','CTCF_10A_L_2','CTCF_10A_L_3','CTCF_10A_L_4','CTCF_10A_L_5','CTCF_10A_L_6','CTCF_10A_L_7','CTCF_10A_L_8','CTCF_10A_L_9','CTCF_10A_L_10','CTCF_10A_R_1','CTCF_10A_R_2','CTCF_10A_R_3','CTCF_10A_R_4','CTCF_10A_R_5','CTCF_10A_R_6','CTCF_10A_R_7','CTCF_10A_R_8','CTCF_10A_R_9','CTCF_10A_R_10']
# chr conversion
chr = df.loc[:, ['chr']]
new_chr = df['chr'].str[-1]
new_chr = new_chr.replace({'X': '23'})
new_chr = pd.to_numeric(new_chr)

# # Add methylation data
new_df['me'] = df[me_cols].values.tolist()
new_df['me'] = new_df['me'].apply(lambda x: np.array(x, dtype=float))
new_df['me'] = new_df['me'].apply(lambda x: np.mean(x))

# # Add histone data
new_df['hist'] = df[hist_cols].values.tolist()
new_df['hist'] = new_df['hist'].apply(lambda x: np.array(x, dtype=float))
new_df['hist'] = new_df['hist'].apply(lambda x: np.mean(x))

# # Add ctcf data
new_df['ctcf'] = df[ctcf_cols].values.tolist()
new_df['ctcf'] = new_df['ctcf'].apply(lambda x: np.array(x, dtype=float))
new_df['ctcf'] = new_df['ctcf'].apply(lambda x: np.mean(x))

# # Add chr
new_df['chr'] = new_chr

# # Add tf data
tf = (df.filter(regex='tfbs_').columns)
new_df['tf'] = df[tf].values.tolist()
new_df['tf'] = new_df['tf'].apply(lambda x: np.array(x, dtype=float))
new_df['tf'] = new_df['tf'].apply(lambda x: np.mean(x))


X = df.loc[:,['me']]

y = df.loc[:, "bds"]
X = np.nan_to_num(X, 0)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

# #print(np.array(X_train))

# #print(np.asarray(np.array(X_train)))


# Running the custom algorithm 
custom = PytorchBasedGenericGradientBoost("classifier", NUM_CLASSIFIERS, MAX_DEPTH, GRADIENT_BOOST_LEARNING_RATE=GRADIENT_BOOST_LEARNING_RATE, MINIMIZER_LEARNING_RATE=MINIMIZER_LEARNING_RATE, MINIMIZER_TRAINING_EPOCHS=MINIMIZER_TRAINING_EPOCHS)
custom.fit(X_train, y_train)
predictions_train = custom.predict(X_train)
predictions_test = custom.predict(X_test)
print("Custom Implementation : Accuracy score for training data : {}".format(accuracy_score(np.round(predictions_train), y_train)))
print("Custom Implementation : Accuracy score for testing data : {}".format(accuracy_score(np.round(predictions_test), y_test)))
torch.save(custom,'./model.pth')

# Running the vanilla sklearn algorithm
# classifier = GradientBoostingClassifier(n_estimators=NUM_CLASSIFIERS, max_depth=MAX_DEPTH)
# classifier.fit(np.array(X_train), np.array(y_train))
# y_pred_gb = custom.predict(X_test)
# print("Vanilla Implementation : Accuracy score for training data : {}".format(accuracy_score(classifier.predict(X_train), y_train)))
# print("Vanilla Implementation : Accuracy score for testing data : {}".format(accuracy_score(classifier.predict(X_test), y_test)))

mae_gb = mean_absolute_error(y_test, y_pred_gb)

mse_gb = mean_squared_error(y_test, y_pred_gb)

importance_gb = custom.feature_importances_
columns = X_train.columns

print(importance_gb)
gbGraph = pd.Series(importance_gb, columns)

print(gbGraph)