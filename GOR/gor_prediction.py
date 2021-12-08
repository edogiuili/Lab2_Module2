from sys import argv
import numpy as np


def read_profile():
   file = open(argv[1],"r")
   submatrix_list = ["H","E","-"]
   row_list = []
   #create the list containing the 51 row names e.i: E|-7
   for subm in submatrix_list:
      for i in range(-8,9):
         row_list.append(subm+"|"+str(i))
   #open the final matrix obtained with gor-prediction and tranform in a numpy matrix
   list_lines_gormatrix = []
   matrix = open("log_matrix.txt","r")
   for line in matrix:
      list_lines_gormatrix.append([float(i) for i in line.split(" ")])

   gor_matrix = np.matrix(list_lines_gormatrix)
   #iterate over the file containing the seq ids
   for seq in file:
      #open a file where write the predicted ss and write the id of the prot
      new_file = open("train_profiles/predicted_dssp/"+seq.rstrip()+".predicted.dssp","w")
      new_file.write(">"+seq+"\n")
      dssp_file = open("train_profiles/"+seq.rstrip()+".dssp","r").readlines()
      list_lines_profile = []
      #call the function open_profile
      profile = open_profile(seq)
      #build a matrix starting from a profile
      for line in profile:
         list_lines_profile.append([float(i.strip()) for i in line.split(" ")])
      profile = np.matrix(list_lines_profile)
      rows,cols = profile.shape
      #iterate over the rows of the profile
      for i in range(rows):
         #initialize the probabilities for each ss
         P_H_Rs = 0
         P_E_Rs = 0
         P_C_Rs = 0
         #iterate over the window
         for j in range(i-8,i+9):
            if j >= 0 and j<rows:
               #iterate over the residues of each position
               for k in range(20):
                  #define the freq of each residue
                  P_R = profile[j,k]
                  log_ratio_H = gor_matrix[row_list.index("H|"+str(j-i)),k]
                  log_ratio_E = gor_matrix[row_list.index("E|"+str(j-i)),k]
                  log_ratio_C = gor_matrix[row_list.index("-|"+str(j-i)),k]
                  #sum the joint probs for each position of the window
                  P_H_Rs += P_R*log_ratio_H
                  P_E_Rs += P_R*log_ratio_E
                  P_C_Rs += P_R*log_ratio_C

         max_prob = max(P_H_Rs,P_C_Rs,P_E_Rs) #use the argmax numpy function
         if max_prob == P_H_Rs:
            new_file.write("H")
         elif max_prob == P_E_Rs:
            new_file.write("E")
         else:
            new_file.write("-")
      #close the file with the predicted ss
      new_file.close()

def open_profile(seq):
   #open the profile
   profile = open("train_profiles/"+seq.rstrip()+".profile.txt","r")
   return profile

read_profile()
