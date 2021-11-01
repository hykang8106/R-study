# to write python version of "meta_test_joon_win10.R"
# test "checkLDPlink" function

'''
fwrite(dat, file = "before_checkLDPlink.csv")
cat("######## before checkLDPlink call saved", sprintf('%s', Sys.time()), "\n")

dat.ld <- checkLDPlink(dat$snpname %>% unique %>% setdiff(c(NA, '')))

# checkLDPlink
checkLDPlink <- function(SNPList) {
  inputFile <- sample(0:9, size = 20, replace = TRUE) %>% paste(collapse = '')
  on.exit(unlink(sprintf('reference/Plink/%s', inputFile)), add = TRUE)
  data.frame(snp = SNPList) %>%
    fwrite(sprintf('reference/Plink/%s', inputFile), quote = FALSE, col.names = FALSE)
  outputFile <- sample(0:9, size = 20, replace = TRUE) %>% paste(collapse = '')
  on.exit(unlink(sprintf('reference/Plink/%s.*', outputFile)), add = TRUE)

  ##### khy
  # must use plink 1.9, DO NOT use plink 1.07
  # for plink 1.9 options,
  # see "https://www.cog-genomics.org/plink/1.9/index"
  # to handle nicely very long format in "sprintf", 
  # see "https://stackoverflow.com/questions/28018473/how-to-avoid-linebreak-in-rs-sprintfvery-very-long-string-with-line-break"
  system(sprintf("C:/Users/hykang/Downloads/plink_win64_20210606/plink --bfile C:/Users/hykang/Desktop/R_study/reference/Plink/ref --memory 1536 --r2 --ld-window-kb 1000 --ld-window 5000 --ld-window-r2 0.7 --ld-snp-list C:/Users/hykang/Desktop/R_study/reference/Plink/%s --out C:/Users/hykang/Desktop/R_study/reference/Plink/%s", 
    inputFile, outputFile), 
    intern = TRUE)
  return(tryCatch(fread(sprintf('reference/Plink/%s.ld', outputFile), encoding = "UTF-8")[, .(SNP_A, SNP_B)], 
                  error = function(e) data.table(SNP_A = character(), SNP_B = character())
                 )
        )
}
'''

import os, sys, random, string, glob
import pandas as pd

from file_db_define import *

test_input_file = "C:/Users/hykang/Desktop/R_study/before_checkLDPlink.csv"

# plink_command_path = 'C:/Users/hykang/Downloads/plink_win64_20210606/plink'
# plink_command_option = \
#  '--memory 1536 --r2 --ld-window-kb 1000 --ld-window 5000 --ld-window-r2 0.7'

# input_file_dir = 'C:/Users/hykang/Desktop/R_study/reference/Plink'
# output_file_dir = 'C:/Users/hykang/Desktop/R_study/reference/Plink'

# plink_ref_path_base = 'C:/Users/hykang/Desktop/R_study/reference/Plink/ref'

def check_LD_Plink()
# generate temporary filename: random 20 digit string
tmp_input_file = ''.join(random.choices(string.digits, k=20))
tmp_output_file = ''.join(random.choices(string.digits, k=20))

input_file_path = '{}/{}'.format(input_file_dir, tmp_input_file)
output_file_path = '{}/{}'.format(output_file_dir, tmp_output_file)

system_string = '{} --bfile {} {} --ld-snp-list {} --out {}'.\
        format(plink_command_path, plink_ref_path_base, \
            plink_command_option, input_file_path, output_file_path)
# print(system_string)

# extract "snpname" column, and make unique
test_input = pd.read_csv(test_input_file, encoding='UTF-8')["snpname"].unique()

# above "test_input" returned data type is "numpy.ndarray", 
# so convert to dataframe, and remove NA
test_input = pd.DataFrame(test_input).dropna()

# write "snpname" dataframe into "input_file_path"
# this is used for input of "plink" command
test_input.to_csv(input_file_path, index=False, header=False)

# execute "plink" command, "plink" make output file 
os.system(system_string)

# read "plink" output file whose name is output_file_path + ".ld"
dat_ld = pd.read_csv(output_file_path + ".ld", encoding='UTF-8', sep='\s+', engine='python')

# extract "SNP_A", "SNP_B" column
dat_ld = dat_ld[['SNP_A', 'SNP_B']]

# remove input file used for "plink" command
try:
  os.remove(input_file_path)
except:
  print("Error while deleting file : ", input_file_path)
# remove output file used for "plink" command
# one more file is related with output

# Get a list of all the file paths using wild character('.*')
fileList = glob.glob(output_file_path + '.*')
# Iterate over the list of filepaths & remove each file.
for filePath in fileList:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)

if __name__ == "__main__":


