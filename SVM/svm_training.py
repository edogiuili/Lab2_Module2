import numpy as np
from sys import argv
from sklearn import svm
import pickle, gzip

def create_matrices(file,directory):
   feature_list = []
   class_list = []
   seq_file = open(file,"r")
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

def define_model(directory,C_list,Y_list,train_set):
   for C in C_list:
      for gamma in Y_list:
         for i in range(len(train_set)):
            train_set = train_cv[i]
            train_class_list,train_feature_list = create_matrices(train_set,directory)
            #create an SVC model
            mySVC = svm.SVC(C=C,kernel='rbf',gamma=gamma)
            #train the SVC model
            print("Fitting the model with C="+str(C)+" and Gamma="+str(gamma)+" in "+train_set)
            mySVC.fit(train_feature_list,train_class_list)
            #save the model
            print("Saving the model...")
            pickle.dump(mySVC,gzip.open(train_set+"_"+"C"+str(C)+"_"+"y"+str(gamma)+".pkl.gz","w"))

directory = argv[1]
C_list = [2,4]
Y_list = [0.5,2]
train_set = [argv[2]]
define_model(directory,C_list,Y_list,train_set)
