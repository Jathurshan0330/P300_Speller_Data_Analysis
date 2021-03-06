# -*- coding: utf-8 -*-
"""P300_Speller.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1r7vk0oO9e72ZjvIGkTP_UfG6hFk0LgOS
"""

!pip install mne

from mne.decoding import CSP
from sklearn.model_selection import GridSearchCV

import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
import scipy.io

from google.colab import drive
drive.mount('/content/drive')

cd '/content/drive/MyDrive/P300_Speller_data/first_experiment/600ms_sw'

pos_data = scipy.io.loadmat('pos_data_S5_600ms_Spectral_whitened.mat')['pos_data']
print(pos_data.shape)
neg_data = scipy.io.loadmat('neg_data_S5_600ms_Spectral_whitened.mat')['neg_data']
print(neg_data.shape)

"""**Training Data**"""

data = np.append(pos_data ,neg_data ,axis = 2) 
data = np.moveaxis(data,-1,0)
data = np.moveaxis(data,-1,1)
data_ = data.reshape(data.shape[0],-1)
print(data_.shape)

labels = np.append(np.ones((pos_data.shape[2],)),np.zeros((neg_data.shape[2],)),axis=0)
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=1)
print(labels.shape)

"""**Classifier**"""

#csp = CSP(n_components=8, reg=None, log=True)

#A = csp.fit_transform(data,labels)

#print(A.shape)

clf = LDA()
scores = cross_val_score(clf, data_, labels, scoring='accuracy', cv=cv, n_jobs=-1)

'''
grid = dict()
grid['solver'] = ['svd', 'lsqr', 'eigen']
# define search
search = GridSearchCV(clf, grid, scoring='accuracy', cv=cv, n_jobs=-1)
# perform the search
results = search.fit(data_, labels)
# summarize
print('Mean Accuracy: %.3f' % results.best_score_)
print('Config: %s' % results.best_params_)
'''

print(np.mean(scores))
print(np.std(scores))

scores

