from sys import argv
import numpy as np
import gzip,pickle

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


def define_model(directory,C_list,Y_list,train_set,test_set):
   performances_dict = {}
   for C in C_list:
      for gamma in Y_list:
         MCC_same_param = 0
         
         test_class_list,test_feature_list = create_matrices(test_set,directory)
         model_file = train_set+"_"+"C"+str(C)+"_"+"y"+str(gamma)+".pkl.gz"
         print(model_file)

         #predict the ss
         print("Predicting the secondary structures in "+test_set)
         mySVC = pickle.load(gzip.open(model_file,"r"))
         ss_pred = mySVC.predict(test_feature_list)

         #compute the performance of this pair of parameters
         print("Computing the performances in "+model_file)
         MCC_avg,MCC_H,MCC_E,MCC_C,Q3,PPV,recall = performance(ss_pred,test_class_list)
         print("MCC_H: ",MCC_H)
         print("MCC_E: ",MCC_E)
         print("MCC_C: ",MCC_C)
         print("Q3: ",Q3)
         print("PPV: ",PPV)
         print("TPR: ",recall)
         MCC_same_param += MCC_avg

      mcc = "MCC_"+str(C)+"_"+str(gamma)
      avg_MCC_same_param = MCC_same_param/len(train_cv)
      print(mcc+":"+str(avg_MCC_same_param)+"\n")
      performances_dict[mcc] = avg_MCC_same_param

   max_key = max(performances_dict, key=performances_dict.get)
   print(max_key)
   print(performances_dict[max_key])

def performance(ss_pred,test_class_list):
   actual_dssp = test_class_list
   predicted_dssp = list(ss_pred)
   class_matrix = np.zeros(shape=(3,3))

   for i in range(len(predicted_dssp)):
      row = actual_dssp[i]-1
      col = predicted_dssp[i]-1
      class_matrix[row,col] += 1

   tp_H = class_matrix[0,0]
   tn_H = class_matrix[1,1]+class_matrix[1,2]+class_matrix[2,1]+class_matrix[2,2]
   fn_H = class_matrix[0,1]+class_matrix[0,2]
   fp_H = class_matrix[1,0]+class_matrix[2,0]

   tp_E = class_matrix[1,1]
   tn_E = class_matrix[0,0]+class_matrix[0,2]+class_matrix[2,0]+class_matrix[2,2]
   fn_E = class_matrix[1,0]+class_matrix[1,2]
   fp_E = class_matrix[0,1]+class_matrix[2,1]

   tp_C = class_matrix[2,2]
   tn_C = class_matrix[0,0]+class_matrix[1,0]+class_matrix[0,1]+class_matrix[1,1]
   fn_C = class_matrix[2,0]+class_matrix[2,1]
   fp_C = class_matrix[0,2]+class_matrix[1,2]

   Q3 = (tp_H+tp_E+tp_C)/len(predicted_dssp)
   mcc_H = ((tp_H*tn_H)-(fp_H*fn_H))/np.sqrt((tp_H+fn_H)*(tp_H+fp_H)*(tn_H+fp_H)*(tn_H+fn_H))
   mcc_E = ((tp_E*tn_E)-(fp_E*fn_E))/np.sqrt((tp_E+fn_E)*(tp_E+fp_E)*(tn_E+fp_E)*(tn_E+fn_E))
   mcc_C = ((tp_C*tn_C)-(fp_C*fn_C))/np.sqrt((tp_C+fn_C)*(tp_C+fp_C)*(tn_C+fp_C)*(tn_C+fn_C))
   mcc_avg = (mcc_H+mcc_E+mcc_C)/3
   recall_H = tp_H/(tp_H+fn_H)
   precision_H = tp_H/(tp_H+fp_H)
   recall_E = tp_E/(tp_E+fn_E)
   precision_E = tp_E/(tp_E+fp_E)
   recall_C = tp_C/(tp_C+fn_C)
   precision_C = tp_C/(tp_C+fp_C)
   precision = (precision_H+precision_E+precision_C)/3
   recall = (recall_H+recall_E+recall_C)/3
   return mcc_avg,mcc_H,mcc_E,mcc_C,Q3,precision,recall

directory = argv[1]
C_list = [2,4]
Y_list = [0.5,2]
#specify the dataset for training
train_set = argv[2]
#specify the dataset for testing
test_set = argv[3]
define_model(directory,C_list,Y_list,train_set,test_set)


