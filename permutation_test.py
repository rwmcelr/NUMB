from random import sample, shuffle
import numpy as np
import pandas as pd
from optparse import OptionParser
import os, sys, time

def import_data(location):
  os.chdir(location)
  dat = pd.read_csv("dat.csv", index_col=0)
  numHigh = dat['UV_sig'].value_counts()['High']
  numLow = dat['UV_sig'].value_counts()['Low']
  return dat.iloc[:,:-1].to_numpy(), numHigh, numLow, dat.iloc[:,:-1].columns

def find_presences(input_matrix):
  num_rows, num_cols = input_matrix.shape
  hp = []
  iters = num_rows if num_cols >= num_rows else num_cols
  input_matrix_b = input_matrix if num_cols >= num_rows else np.transpose(input_matrix)
  for r in range(iters):
    hp.append(list(np.where(input_matrix_b[r] ==1)[0]))
  return hp

def curve_ball(input_matrix, r_hp, num_iterations=-1):
  num_rows, num_cols = input_matrix.shape
  l = range(len(r_hp))
  num_iters = 5*min(num_rows, num_cols) if num_iterations == -1 else num_iterations
  for rep in range(num_iters):
    AB = sample(l, 2)
    a = AB[0]
    b = AB[1]
    ab = set(r_hp[a])&set(r_hp[b])
    l_ab=len(ab)
    l_a=len(r_hp[a])
    l_b=len(r_hp[b])
    if l_ab not in [l_a,l_b]:
      tot=list(set(r_hp[a]+r_hp[b])-ab)
      ab=list(ab)
      shuffle(tot)
      L=l_a-l_ab
      r_hp[a] = ab+tot[:L]
      r_hp[b] = ab+tot[L:]
  out_mat = np.zeros(input_matrix.shape, dtype='int8') if num_cols >= num_rows else np.zeros(input_matrix.T.shape, dtype='int8')
  for r in range(min(num_rows, num_cols)):
    out_mat[r, r_hp[r]] = 1
  result = out_mat if num_cols >= num_rows else out_mat.T
  return result

def permTest(options):
  M,h,l,colNames = import_data(options.data_location)
  r_hp = find_presences(M)
  M_dif = M[0:h,:].sum(axis=0)-M[h:l+h,:].sum(axis=0)
  results = pd.DataFrame(0, index=colNames, columns=["p"])
  start = time.time()
  for i in range(options.num_perms):
    RM = curve_ball(M, r_hp)
    RM_dif = RM[0:h,:].sum(axis=0)-RM[h:l+h,:].sum(axis=0)
    for j in range(len(M_dif)):
      if RM_dif[j] >= M_dif[j]:
        results.iloc[[j]] += 1
    if i % 100 == 0:
      print("Generated {0} matrices in {1} seconds, {2} permutations remaining.".format(i, time.time()-start, options.num_perms-i))
  (results/options.num_perms).to_csv(str(options.num_perms)+'_Permutations_Raw.csv', index=True)

def main():
  parser = OptionParser()
  parser.add_option("-d", "--data_location", dest="data_location", help="Matrix input file")
  parser.add_option("-n", "--num", dest="num_perms", default=1000, type="int", help="Number of Permutations to run")

  (options, args) = parser.parse_args()

  # Run permutation test
  permTest(options)

if __name__ == '__main__':
  sys.exit(main())
