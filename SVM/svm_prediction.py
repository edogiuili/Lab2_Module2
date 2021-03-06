from sys import argv
import numpy as np
import gzip,pickle
import sklearn.metrics 


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


def define_model(directory,C_list,Y_list,train_cv,test_cv):
   performances_dict = {}
   perf_file = open("performances.txt","w")
   for C in C_list:
      for gamma in Y_list:
         MCC_same_param = 0
         for i in range(len(train_cv)):
            train_set = train_cv[i]
            test_set = test_cv[i]
            test_class_list,test_feature_list = create_matrices(test_set,directory)
            model_file = train_set+"_"+"C"+str(C)+"_"+"y"+str(gamma)+".pkl.gz"

            #predict the ss
            print("Predicting the secondary structures in "+test_set)
            mySVC = pickle.load(gzip.open(model_file,"r"))
            ss_pred = mySVC.predict(test_feature_list)

            #compute the performance of this pair of parameters
            print("Computing the performances in "+model_file)
            Q3,MCC_avg,MCC_H,MCC_E,MCC_C,acc_H,acc_E,acc_C,precision_H,precision_E,precision_C,recall_H,recall_E,recall_C = performance(ss_pred,test_class_list)
            perf_file.write("Performances of "+model_file+"\n")
            perf_file.write("MCC_H:"+str(MCC_H)+" ")
            perf_file.write("MCC_E:"+str(MCC_E)+" ")
            perf_file.write("MCC_C:"+str(MCC_C)+" ")
            perf_file.write("Q3:"+str(Q3)+" ")
            perf_file.write("Acc_H:"+str(acc_H)+" ")
            perf_file.write("Acc_E:"+str(acc_E)+" ")
            perf_file.write("Acc_C:"+str(acc_C)+" ")
            perf_file.write("PPV_H:"+str(precision_H)+" ")
            perf_file.write("PPV_E:"+str(precision_E)+" ")
            perf_file.write("PPV_C:"+str(precision_C)+" ")
            perf_file.write("TPR_H:"+str(recall_H)+" ")
            perf_file.write("TPR_E:"+str(recall_E)+" ")
            perf_file.write("TPR_C:"+str(recall_C)+"\n")
            MCC_same_param += MCC_avg

         mcc_model_name = "MCC_"+str(C)+"_"+str(gamma)
         avg_MCC_same_param = MCC_same_param/len(train_cv)
         performances_dict[mcc_model_name] = avg_MCC_same_param

   max_key = max(performances_dict, key=performances_dict.get)
   print(max_key)
   print(performances_dict[max_key])
   perf_file.close()

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

   acc_H = (tp_H+tn_H)/(tn_H+tp_H+fp_H+fn_H)
   acc_E = (tp_E+tn_E)/(tn_E+tp_E+fp_E+fn_E)
   acc_C = (tp_C+tn_C)/(tn_C+tp_C+fp_C+fn_C)

   Q3 = (tp_H+tp_E+tp_C)/class_matrix.sum()

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
   return Q3,mcc_avg,mcc_H,mcc_E,mcc_C,acc_H,acc_E,acc_C,precision_H,precision_E,precision_C,recall_H,recall_E,recall_C

directory = argv[1]
C_list = [2,4]
Y_list = [0.5,2]
train_set = ["train0","train1","train2","train3","train4"]
test_set = ["test4","test3","test2","test1","test0"]
define_model(directory,C_list,Y_list,train_set,test_set)


