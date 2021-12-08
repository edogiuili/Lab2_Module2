from sys import argv
import numpy as np


def create_matrix():
   residues_list = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
   submatrix_list = ["H","E","-","R"]
   marginal_prob_ss = {"H":0,"E":0,"-":0}
   row_list = []
   #create the list containing the 68 row names e.i: E|-7
   for subm in submatrix_list:
      for i in range(-8,9):
         row_list.append(subm+"|"+str(i))
 
   #initialize the matrix with all zeros and add the residue row and the last column
   matrix = np.zeros(shape=(68,20))
   return(row_list,matrix,marginal_prob_ss)

def read_profile(row_list,matrix,marginal_prob_ss,file,directory):
   #initialize the counter for the final normalization
   total_length = 0
   #iterate over the seqs contained in the match_seqs file
   for seq in file:
      #define the profile and the dssp sequence calling iter_files function
      try:
         profile, dssp_sequence = iter_files(directory,seq)
         profile = profile.readlines()
         dssp_sequence = dssp_sequence.readlines()[1].rstrip()
         #iterate over the dssp sequence
         for i in range(len(dssp_sequence)):
            #let be the secondary structure of the element the ss for all the window
            ss = dssp_sequence[i]
            #iterate over the window
            for j in range(i-8,i+9):
               #take only positions in the sequence
               if j >= 0 and j < len(dssp_sequence):
                  row = row_list.index(ss+"|"+str(j-i))
                  marginal_rowres = row_list.index("R|"+str(j-i))
                  #iterate over the profile residues (column index)
                  for k in range(20):
                     #add the profile scores to the gor matrix
                     matrix[row,k] += float(profile[j].split()[k])
                     matrix[marginal_rowres,k] += float(profile[j].split()[k])
            #add 1 to the marginal probs
            marginal_prob_ss[ss] += 1
         #add the length of the sequence to the total length
         total_length += len(dssp_sequence)
      except Exception as message:
         print(message)
   #normalize the matrix for the total length
   matrix = matrix/total_length
   #normalize the marginal probs of the ss
   for key in marginal_prob_ss.keys():
      marginal_prob_ss[key]=marginal_prob_ss[key]/total_length
   file.close()
   return matrix, marginal_prob_ss,total_length


def iter_files(directory,seq):
   profile_seq = open(directory+seq.rstrip()+".profile.txt","r")
   dssp_seq = open(directory+seq.rstrip()+".dssp","r")
   return profile_seq,dssp_seq

def log_matrix(gor_matrix):
   #define the dimensions of the matrix
   rows,cols=gor_matrix.shape
   #iterate over the rows (positions) of the matrix
   for i in range(rows):
      #iterate over the 20 residues
      for j in range(cols):
         #compute the log score for each position of the matrix
         if i >= 0 and i < 17:
            gor_matrix[i,j] = np.log(gor_matrix[i,j]/(gor_matrix[i+51,j]*marginal_prob_ss["H"]))
         if i >= 17 and i < 34:
            gor_matrix[i,j] = np.log(gor_matrix[i,j]/(gor_matrix[i+34,j]*marginal_prob_ss["E"]))
         if i >= 34 and i < 51:
            gor_matrix[i,j] = np.log(gor_matrix[i,j]/(gor_matrix[i+17,j]*marginal_prob_ss["-"]))
   return(gor_matrix)

def write_matrix(log_gor_matrix,output_file):
   with open(argv[2]+argv[3],"w") as f:
      gor_matrix_list = gor_matrix[:51].tolist()
      for position in gor_matrix_list:
         line_str = [str(i) for i in position]
         f.write(" ".join(line_str)+"\n")
   f.close()

file = open(argv[1],"r")
directory = argv[2]
output_file = argv[3]

row_list,matrix,marginal_prob_ss = create_matrix()
matrix_final,ss_prob_final,total_length = read_profile(row_list,matrix,marginal_prob_ss,file,directory)
log_gor_matrix = log_matrix(matrix_final)
write_matrix(log_gor_matrix,output_file)
