## Example of training run
*python3 gor_training.py <input_train_seqs> <directory> <output_matrix_file>*
 
The directory must contain all the training sequences profiles
  
## Example of predicting run
  *python3 gor_prediction.py <input_test_seqs> <directory> <log_matrix_file>
  The directory must contain all the test set sequences profiles

## Example of performance run
*python3 gor_performance.py <input_test_seqs> <directory>
  The directory must contain all the test set dssp sequences (predicted in the previous run and expected)
