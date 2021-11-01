

'''
  (function(dat) {
    # load reference bim from file
    ##### khy
    # "ref.bim" file have data.table
    # set data.table("bim") column name to "CHR", "snpname", "POS", "REF", "ALT"
    fread('reference/Plink/ref.bim', encoding = "UTF-8", header = FALSE)[, .(CHR = V1, snpname = V2, POS = V4, REF = V5, ALT = V6)] %>%  
      (function(bim) {
        # create "SNP" column, and input rs number extracted from "snpname" column
        # for ":=" operator, 
        # see "https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reference-semantics.html"
        bim[, SNP := gsub('^(rs\\d+)(\\D.*)?$', '\\1', snpname)] 
        # set key "SNP" column
        setkey(bim, SNP)
        # "bim[dat]" is DT syntax, which is same as "merge(bim, dat)"
        # see DT syntax:
        # [ref] https://rstudio-pubs-static.s3.amazonaws.com/52230_5ae0d25125b544caab32f75f0360e775.html
        afterbim <- bim[dat]
        fwrite(afterbim, file = "afterbim.csv")
        #fwrite(afterbim, file = "afterbim.rdata")
        #### khy
        # to seee "afterbim" in rterm, use "View(fread("afterbim.rdata"))"
        cat("##### after bim saved", sprintf('%s', Sys.time()), "\n")
        afterbim
      }) 
  }) %>%
'''

'''
# read "g1" from temporary file
g1 = pd.read_csv(before_bim_file, encoding="UTF-8")
print(g1.shape)
print(g1.head())


#### ref_bim_file not have column name, separation character = "tab", column name prefix = 'V'
bim = pd.read_csv(ref_bim_file, encoding="UTF-8", header=None, sep='\t', prefix='V')
# print(bim.shape)
# print(bim.head())

# select column from "bim" dataframe: R = one-indexing, python = zero-indexing
bim = bim[['V0', 'V1', 'V3', 'V4', 'V5']]

# rename column name
bim = bim.rename(columns={'V0':'CHR', 'V1':'snpname', 'V3':'POS', 'V4':'REF', 'V5':'ALT'})
print(bim.shape)
print(bim.head())

# add 'SMP' column, and fill with replaced data from 'snpname' column: "rs10001545:C:A" => "rs10001545"
# [ref] https://stackoverflow.com/questions/37425019/gsub-only-part-of-pattern
# '\D' meaning: matches a non-digit
bim['SNP'] = [re.sub('^(rs\d+)(\D.*)?$', '\\1', s) for s in bim['snpname']]
print(bim.head())

# bim.to_csv('before_merge_bim_py.csv', encoding='utf-8', index=False)

# merge "bim" with "g1": on = 'SNP', how = 'right'('right' is right?)
# how option = ‘left’, ‘right’, ‘outer’, ‘inner’
# [ref] https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html
# [ref] https://stackoverflow.com/questions/53645882/pandas-merging-101/53645883#53645883
bim = bim.merge(g1, on='SNP', how='right')
#### check: "merge" convert 'CHR', 'POS' colmun from integer to float, why?

bim.to_csv('afterbim_py.csv', encoding='utf-8', index=False)
'''

'''
  # if "CHR"(Chromosome) column is empty, remove that row
  filter(!is.na(CHR)) %>% 
  # if "eff_allele" is not same as "REF" or "ALT", remove that row
  filter(REF == eff_allele | ALT == eff_allele) %>% 
  #filter(eff_allele == REF | eff_allele == ALT) %>% 
  # remove row having duplicated "Bio.Serial"(it is due to multi allele?)
  filter(!duplicated(Bio.Serial)) %>%
  # create "ref_allele" column, and
  # if "eff_allele" is same as "ALT", input "REF", otherwise "ALT" column
  mutate(ref_allele = ifelse(eff_allele == ALT, REF, ALT)) %>% 
  # create "aa" column, and input "ref_allele""ref_allele"
  mutate(aa = sprintf('%s%s', ref_allele, ref_allele)) %>% 
  # create "ab" column, and input "ref_allele""eff_allele" or "eff_allele""ref_allele"
  mutate(ab = ifelse(eff_allele > ref_allele, 
                     sprintf('%s%s', ref_allele, eff_allele),
                     sprintf('%s%s', eff_allele, ref_allele))) %>%
  # create "bb" column, and input "eff_allele""eff_allele" column
  mutate(bb = sprintf('%s%s', eff_allele, eff_allele)) %>% 
  ##### khy ???
  # what below "ungroup" do?
  ungroup %>%
  data.table %>%
'''

'''
# keep row only where 'CHR' column is not null
bim = bim[bim['CHR'].notnull()].reset_index(drop=True)

# "merge" converted 'CHR', 'POS' colmun from integer to float, so restote to integer
bim = bim.astype({'CHR':'int64', 'POS':'int64'})

# print(bim.head())

# bim.to_csv('after_bim_py.csv', encoding='utf-8', index=False)

# keep rows only where 'REF' column is same as 'eff_allele' column or 'ALT' column is same as 'eff_allele'
bim_valid = (bim['REF'] == bim['eff_allele']) | (bim['ALT'] == bim['eff_allele'])
bim = bim[bim_valid].reset_index(drop=True)

# remove rows where 'Bio.Serial' column is duplicated
bim = bim.drop_duplicates(subset=['Bio.Serial']).reset_index(drop=True)

bim['ref_allele'] = np.where(bim['eff_allele'] == bim['ALT'], bim['REF'], bim['ALT'])

bim['aa'] = bim['ref_allele'] + bim['ref_allele']

bim['ab'] = np.where(bim['eff_allele'] > bim['ref_allele'], \
  bim['ref_allele'] + bim['eff_allele'], bim['eff_allele'] + bim['ref_allele'])

bim['bb'] = bim['eff_allele'] + bim['eff_allele']
'''


'''
  ###### khy
  # %dopar% -> %do%
  mutate(gene = foreach(snp = SNP, .combine = c, .multicombine = TRUE) %do% {
    con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.gene.grch37.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    res <- dbGetQuery(con, sprintf('select GENE from gene where SNP="%s"', snp))
    if (nrow(res) == 0)
      return('intergenic')
    else
      return(res$GENE[1])
  }) %>%
  mutate(LINK = ifelse(is.na(LINK), '', LINK)) %>%
  ###### khy
  # %dopar% -> %do%
  mutate(TIMES = foreach(row = 1:n(), .combine = c) %do% {
    return(!is.na(BETA[row]) &
             (!code[row] %in% code.continuous.list |
                (sum(code.continuous$code == code[row] &
                      code.continuous$LINK == LINK[row] &
                      !(is.na(code.continuous$mean) | code.continuous$mean == 0)
                    ) > 0
                )
             )
          )
  }) %>%
  ################### khy
  ######## for debug
(function(dat) {
  fwrite(dat, file = "mutate_TIMES.csv")
  #fwrite(dat, file = "mutate_TIMES.rdata")
  cat("##### mutate TIMES saved", sprintf('%s', Sys.time()), "\n")
  dat
}) %>%
'''

import pandas as pd
import numpy as np
# re: regular expression
import re, sys

# from cytogenic_band import *
from query_dbsnp import *

from mutate_foreach import *

from file_db_define import *

def process_bim(bim, continuous_file, code_continuous):

  # G1_TIMES_file = "C:/Users/hykang/Desktop/R_study/G1_TIMES_py.csv"
  # before_merge_g1_samples_file = "C:/Users/hykang/Desktop/R_study/before_merge_g1_samples_py.csv"

  code_continuous_list = pd.read_csv(continuous_file, encoding="UTF-8")['code'].unique()
  # code_continuous_list = pd.read_csv(continuous_file, encoding="UTF-8").code.unique()
  # "unique()" return numpy array
  # below is comment out because conversion to dataframe give wrong result:
  # when "code" is 'CGH', "code not in code_continuous_list" in "bim.util.py" give "True" (wrong result)
  # code_continuous_list = pd.DataFrame(code_continuous_list)
  print(code_continuous_list.shape)
  # numpy array not support "head()""
  # print(code_continuous_list.head())

  for i in range(len(bim)):

    bim.loc[i, 'gene'] = findGeneDbsnp(bim['SNP'][i])

  # need below?
  bim['LINK'] = np.where(bim['LINK'].isnull(), '', bim['LINK'])

  for i in range(len(bim)):

    bim.loc[i, 'TIMES'] = mutate_TIMES_column(code_continuous_list, code_continuous, \
      bim['BETA'][i], bim['code'][i], bim['LINK'][i])

  bim.to_csv(mutate_TIMES_file, encoding='utf-8', index=False)
  # bim.to_csv('mutate_TIMES_py.csv', encoding='utf-8', index=False)

  # sys.exit()
  '''
    mutate(BETA = foreach(row = 1:n(), .combine = c) %do% {
      if (sum(TIMES[code == code[row]]) == 0)
        return(0.001 * BETA[row] / abs(BETA[row]))
      else if (!code[row] %in% code.continuous.list)
        return(BETA[row])
      mean <- code.continuous$mean[code.continuous$code == code[row] &
                                    code.continuous$LINK == LINK[row]]
      ##### khy
      # "mean" column have non-digit web link
      # this result from corrupted "code.continuous.RData"
      # to see it, use "View(fread("code.continuous.RData"))" or "load('code.continuous.RData')"
      # "mean" column in code.continuous.RData" have many empty!
      mean <- mean[check.numeric(mean)]
      # convert to numneric
      mean <- as.numeric(mean)
      mean <- mean %>% setdiff(NA)
      if (length(mean) == 0) {
        #### khy, for debug
        cat("$$$$ if (length(mean) == 0)\n")
        return(0.001 * BETA[row] / abs(BETA[row]))
      }
      mean <- mean[1]
      if (mean == 0) {
        #### khy, for debug
        cat("$$$$ if (mean == 0)\n")
        return(0.001 * BETA[row] / abs(BETA[row]))
      }
      else if (mean < 0) {
        #### khy, for debug
        cat("$$$$ else if (mean < 0)\n")
        return(BETA[row])
      }
      else {
        #### khy, for debug
        # result: 13 case printed
        #cat("$$$$ else: return(log((BETA[row] + mean) / mean))\n")
        return(log((BETA[row] + mean) / mean))
      }
        ##### khy
        # "log" must have positive value
        #if ((BETA[row] + mean) / mean > 0)
        #  return(log((BETA[row] + mean) / mean))
        #else
        #  return(NA)
    }) %>%
    filter(!is.na(BETA)) %>%
  '''

  for i in range(len(bim)):
    # print(i)
    bim.loc[i, 'BETA'] = mutate_BETA_column(code_continuous_list, code_continuous, bim, i)

  beta_valid = bim['BETA'].notnull()
  bim = bim[beta_valid].reset_index(drop=True)

  '''
  (function(dat) {
    dat %>%
      group_by(code) %>%
      #### khy
      # below "na.rm = TRUE" mean NaN removed
      # used to count TRUE values in logical vector
      summarize(TIMES = sum(TIMES, na.rm = TRUE)) %>%
      filter(TIMES == 0) %>%
      dplyr::select(code) %>%
      fwrite('G1_TIMES.csv')
      cat("##### G1_TIMES saved", sprintf('%s', Sys.time()), "\n")
    dat
  }) %>%
  '''

  write_G1_TIMES(bim, G1_TIMES_file)
  ####### what "G1_TIMES" is for? "G1_TIMES" is empty. right? 
  ####### R code produce non empty "G1_TIMES"
  print("### G1_TIMES saved")

  '''
    mutate(SNP_OR = exp(BETA)) %>%
    (function(dat) {
      dat %>%
        (function(x) {
          ##### khy, debug
          fwrite(x, file = "before_grepl_column_name.csv")
          #fwrite(dat, file = "before_checkLDPlink.rdata")
          cat("######## before grepl column name saved", sprintf('%s', Sys.time()), "\n")
          #### khy
  '''

  bim['SNP_OR'] = np.exp(bim['BETA'])

  bim.to_csv(before_merge_g1_samples_file, encoding='utf-8', index=False)

  # below R code is not translated because it seems to be not necessary

  '''
          # "x %>% names" means to get column names of x object
          # below meaning: 
          # select data in column whose column name is matched with grepl pattern from x
          # "with = FALSE" give all data in column, "with = TRUE" give only column name
          #### khy
          # why below line is needed? there is no column name matched with '\\.(Initial|Replication)$' 
          # my guess: below 'reference/g1_samples.csv' file have column name matched with the pattern
          x[, !grepl('\\.(Initial|Replication)$', x %>% names), with = FALSE]
        }) %>%
  '''

  return bim

if __name__ == "__main__":

  import pandas as pd

  from get_code_continuous_mean import *

  code_continuous = get_code_continuous_mean(g1_continuous_mean_file)

  bim = pd.read_csv(cytogenic_refered_file, encoding='utf-8')

  bim = process_bim(bim, continuous_file, code_continuous)
