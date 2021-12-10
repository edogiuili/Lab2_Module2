import numpy as np
from sys import argv
from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.model_selection import PredefinedSplit
import pickle, gzip

def create_matrices(file,directory):
   feature_list = []
   class_list = []
   seq_file = open(directory+file,"r")
   for seq in seq_file:
      try:
         profile = open(directory+seq.rstrip()+".profile.txt").readlines()
         dssp = open(directory+seq.rstrip()+".dssp").readlines()
         seq_dssp = dssp[1].rstrip()
         #transform the profile in a numpy matrix
         matrix = np.array(profile[0].rstrip().split(" "))
         for line in profile[1:]:
            matrix = np.vstack([matrix,line.rstrip().split(" ")])
         matrix = matrix.astype(np.float64)
         for i in range(len(seq_dssp)):
            position_list = []
            if seq_dssp[i] == "H":
               class_list.append(1)
            elif seq_dssp[i] == "E":
               class_list.append(2)
            elif seq_dssp[i] == "-":
               class_list.append(3)
            for j in range(i-8,i+9):
               if j < 0 or j >= len(seq_dssp):
                  position_list += [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
               else:
                  position_list += list(matrix[j])
            feature_list += [position_list]
      except:
         pass
   return class_list,feature_list

def define_model(directory,train_set,parameters):
   X = []
   Y = []
   for train_split in train_set:
      train_class_list,train_feature_list = create_matrices(train_split,directory)
      X.append(train_feature_list)
      Y.append(train_class_list)
   
   test_fold = [0,1,2]
   cv = PredefinedSplit(test_fold=test_fold)
   #create an SVC model
   model = GridSearchCV(estimator=SVC(),cv=,param_grid=parameters,return_train_score=True)
   #train the SVC model
   model.fit(X,Y)
   #get best accuracy and hyperparms after the grid search
   print(model.best_accuracy_)
directory = argv[1]
parameters = {"C":[2,4],"gamma":[0.5,2],"kernel":["rbf"]}
train_set = ["prova.txt","prova2.txt","prova3.txt"]
define_model(directory,train_set,parameters)
