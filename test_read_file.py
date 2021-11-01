# to write python version of "meta_test_joon_win10.R"
# test read file

import pandas as pd

g1_continuous_mean_file = "C:/Users/hykang/Desktop/R_study/reference/fixed_G1_Continuous_Mean.csv"
plink_ref_freq_file = "C:/Users/hykang/Desktop/R_study/reference/Plink/ref.frq"
continuous_file = "C:/Users/hykang/Desktop/R_study/reference/continuous.csv"
g1_file = "C:/Users/hykang/Desktop/R_study/reference/fixed_g1.csv"
ref_bim_file = "C:/Users/hykang/Desktop/R_study/reference/Plink/ref.bim"
g1_samples_file = "C:/Users/hykang/Desktop/R_study/reference/g1_samples.csv"

## "cytogenic.band.length.RData" file is loaded by "load" function in R:
# "load('reference/cytogenic.band.length.RData')"
# this file may be saved using "save" function in R
# "save" in R may use compression. so in python, "pd.read_csv" not work
# to cure this problem:
# in R, load this file, and write back to file using "fwrite" function (default sep=','):
# "fwrite(cytogenic.band.length, "reference/cytogenic.band.length.csv")"
# then in python, "pd.read_csv" work
# but there may be nicer solution
# try to use "rpy2" package,
# not easy to install in windows 10: pip3 install rpy2 ===> failed
# [ref] https://stackoverflow.com/questions/21288133/loading-rdata-files-into-python
# [ref] https://rpy2.github.io/doc/latest/html/index.html

cytogenic_band_length_file = "C:/Users/hykang/Desktop/R_study/reference/cytogenic.band.length.csv"

g1_continuous_mean = pd.read_csv(g1_continuous_mean_file, encoding="cp949")
print(g1_continuous_mean.shape)
print(g1_continuous_mean.head())

# separation char is one more space, so sep='\s+', and grep expression need 'python' engine
plink_ref_freq = pd.read_csv(plink_ref_freq_file, encoding="UTF-8", sep='\s+', engine='python')
print(plink_ref_freq.shape)
print(plink_ref_freq.head())

continuous = pd.read_csv(continuous_file, encoding="UTF-8")
print(continuous.shape)
print(continuous.head())

g1 = pd.read_csv(g1_file, encoding="cp949")
print(g1.shape)
print(g1.head())

#### ref_bim_file not have column name, separation character = "tab", column name prefix = 'V'
ref_bim = pd.read_csv(ref_bim_file, encoding="UTF-8", header=None, sep='\t', prefix='V')
print(ref_bim.shape)
print(ref_bim.head())

g1_samples = pd.read_csv(g1_samples_file, encoding="UTF-8")
print(g1_samples.shape)
print(g1_samples.head())

cytogenic_band_length = pd.read_csv(cytogenic_band_length_file, encoding="UTF-8")
print(cytogenic_band_length.shape)
print(cytogenic_band_length.head())

