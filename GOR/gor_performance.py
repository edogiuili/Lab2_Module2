from sys import argv
import numpy as np

def compare_files(seq_file,directory):
   #create the confusion matrix 3x3
   class_matrix = np.zeros(shape=(3,3))
   #initialize a dict with the index for each ss
   class_dict = {"H":0,"E":1,"-":2}
   #iterate over the seqs in the file
   for line in seq_file:
      try:
         #fill the confusion matrix for each seq
         predicted_dssp = open(directory+line.rstrip()+".predicted.dssp","r").readlines()
         actual_dssp = open(directory+line.rstrip()+".dssp","r").readlines()
         for i in range(len(predicted_dssp[1].rstrip())):
            row = class_dict[actual_dssp[1][i]]
            col = class_dict[predicted_dssp[1][i]]
            class_matrix[row,col] += 1
      except:
         pass

   #compute the tp,tn,fp,fn for each ss
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

   mcc_H = ((tp_H*tn_H)-(fp_H*fn_H))/np.sqrt((tp_H+fn_H)*(tp_H+fp_H)*(tn_H+fp_H)*(tn_H+fn_H))
   mcc_E = ((tp_E*tn_E)-(fp_E*fn_E))/np.sqrt((tp_E+fn_E)*(tp_E+fp_E)*(tn_E+fp_E)*(tn_E+fn_E))
   mcc_C = ((tp_C*tn_C)-(fp_C*fn_C))/np.sqrt((tp_C+fn_C)*(tp_C+fp_C)*(tn_C+fp_C)*(tn_C+fn_C))
   q3 = (class_matrix[0,0]+class_matrix[1,1]+class_matrix[2,2])/class_matrix.sum()
   recall = (tp_E+tp_C+tp_H)/(tp_H+fn_H+tp_E+fn_E+tp_C+fn_C)
   precision = (tp_E+tp_C+tp_H)/(tp_H+fp_H+tp_E+fp_E+tp_C+fp_C)
   return mcc_H,mcc_E,mcc_C,q3,recall,precision

seq_file = open(argv[1],"r")
directory = argv[2]
MCC_H,MCC_E,MCC_C,Q3,TPR,PPV = compare_files(seq_file,directory)

print("MCC_H : ",round(MCC_H,3))
print("MCC_E : ",round(MCC_E,3))
print("MCC_C : ",round(MCC_C,3))
print("Q3 : ",round(Q3,3))
print("TPR: ",round(TPR,3))
print("PPV: ",round(PPV,3))


