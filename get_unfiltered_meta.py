# to write python version of "meta_test_joon_win10.R"
# get "unfiltered_meta" data

'''
load('reference/cytogenic.band.length.RData')

# process bio team data
#### khy
# 'code.continuous.RData' was created
# see "save(code.continuous, file = 'code.continuous.RData')"
load('code.continuous.RData')
code.continuous.list <- fread('reference/continuous.csv', encoding = "UTF-8")$code %>% unique
'''

import pandas as pd
import numpy as np
# re: regular expression
import re, sys
# import function: "findAltDbsnp", "checkAlleleDbsnp", "findGeneDbsnp", "rev_strand"
from query_dbsnp import *

cytogenic_band_length_file = "C:/Users/hykang/Desktop/R_study/reference/cytogenic.band.length.csv"
code_continous_file = "C:/Users/hykang/Desktop/R_study/code_continuous.csv"
continuous_file = "C:/Users/hykang/Desktop/R_study/reference/continuous.csv"
g1_file = "C:/Users/hykang/Desktop/R_study/reference/fixed_g1.csv"
ref_bim_file = "C:/Users/hykang/Desktop/R_study/reference/Plink/ref.bim"

cytogenic_band_length = pd.read_csv(cytogenic_band_length_file, encoding="UTF-8")
print(cytogenic_band_length.shape)
print(cytogenic_band_length.head())

code_continous = pd.read_csv(code_continous_file, encoding="UTF-8")
print(code_continous.shape)
print(code_continous.head())

code_continuous_list = pd.read_csv(continuous_file, encoding="UTF-8").code.unique()
# "unique()" return numpy array
code_continuous_list = pd.DataFrame(code_continuous_list)
print(code_continuous_list.shape)
print(code_continuous_list.head())

'''
unfiltered_meta <- fread('reference/fixed_g1.csv') %>%
  # Bio.Serial column: Assign serial Number for each SNP
  mutate(Bio.Serial = 1:n()) %>%
  filter(!is.na(code)) %>%
  # Checking history from BS : 1st check is '0'
  filter(`1st check` == 'O') %>%
  
  #### khy
  # what is "RSID %in% c('')"?
  # Substitute RSID col. with SNPS col., if RSID col. is empty and SNPS col. starts with 'rs' and ends with some digits    
  mutate(RSID = ifelse(RSID %in% c('') & grepl('^rs\\d+$', SNPS), SNPS, RSID)) %>%
  # Filter in RSID col. with rs number format
  filter(grepl('rs\\d+', RSID)) %>%
  
  # Rename column as readable in R
  rename(SNP = RSID, 
         sample.EAS = 'Asia(EAS)', 
         sample.OTHER = 'Non-EAS', 
         sample.TOTAL = 'Total sample size') %>%
  # Convert sample.EAS, sample.OTHER, sample.TOTAL col. to numeric value and NAs to 0
  mutate(sample.EAS = ifelse(is.na(sample.EAS), 0, sample.EAS) %>% as.numeric,  
         sample.OTHER = ifelse(is.na(sample.OTHER), 0, sample.OTHER) %>% as.numeric, 
         sample.TOTAL = ifelse(is.na(sample.TOTAL), 0, sample.TOTAL) %>% as.numeric) %>%
  
  ##### khy
  # remove rows whose 'P-VALUE' have string non-convertible to numeric
  # "check.numeric" function is used, which require "varhandle" package
  filter(check.numeric(`P-VALUE`)) %>%

  # P-value to numeric (warning:NAs)
  mutate(`P-VALUE` = as.numeric(`P-VALUE`)) %>%
  
  # Select columns required to calculate
  dplyr::select(code, SNP, sample.EAS, sample.OTHER, sample.TOTAL,
                `OR/RR/BETA`, `OR or BETA`, Effect, `P-VALUE`, PUBMEDID, LINK, `Do/Re`, Bio.Serial) %>%
  
  # Rename column as readable in R
  ##### khy ????
  # 'OR or BETA' = 'OR or BETA', same name was renamed, why?
  rename('OR or BETA' = 'OR or BETA', 'eff_allele' = 'Effect') %>%
  # Unify other form of effect of SNP to the word 'beta'; maintain 'OR'
  #OR/RR/BETA의 정보를 beta 또는 OR으로 변경
  #### khy ????
  # above comment is right? below line show changing to 'BETA' or 'OR/RR/BETA'
  # what is "grepl('(beta|% effect|coef|AUC|delta|mean|estimate|h2|u|z|ln\\(OR\\))'"?
  mutate(`OR/RR/BETA` = ifelse(!is.na(`OR/RR/BETA`) & 
                                grepl('(beta|% effect|coef|AUC|delta|mean|estimate|h2|u|z|ln\\(OR\\))', 
                                  `OR/RR/BETA`, ignore.case = TRUE), 
                          'BETA', `OR/RR/BETA`)) %>% 
  # OR or BETA to numeric (warnings:NAs)

  ##### khy
  filter(check.numeric(`OR or BETA`)) %>%
  mutate(`OR or BETA` = as.numeric(`OR or BETA`)) %>%
  
  # OR to log(OR)
  ##### khy
  # `OR or BETA` have negative value, log of negative make NaN
  filter(`OR or BETA` > 0) %>%
  mutate(BETA = ifelse(!is.na(`OR or BETA`) & grepl('beta', `OR/RR/BETA`, ignore.case = TRUE), 
                        `OR or BETA`, log(`OR or BETA`))) %>%
  
  # publication reference as pubmed link
  # change data in "PUBMEDID" column to integer
  ##### khy
  # "PUBMEDID" have non-digit data
  filter(check.numeric(PUBMEDID)) %>%
  mutate(PUBMEDID = PUBMEDID %>% as.integer) %>%  
  # if "PUBMEDID" is valid, set "LINK" to url name including "PUBMEDID"
  mutate(LINK = ifelse(!is.na(PUBMEDID) &
                         PUBMEDID != '' &
                         PUBMEDID != '0', 
                         sprintf('www.ncbi.nlm.nih.gov/pubmed/%s', PUBMEDID), LINK)) %>% 
  
  # ########################################################################
  # make "DOMINANT" column, and if "Do/Re" column is "Do", set TRUE
  mutate(DOMINANT = grepl('Do', `Do/Re`, ignore.case = TRUE)) %>% 
  # make "RECESSIVE" column, and if "Do/Re" column is "Re", set TRUE 
  mutate(RECESSIVE = grepl('Re', `Do/Re`, ignore.case = TRUE)) %>% 
  # convert "eff_allele" column to upper character 
  mutate(eff_allele = toupper(eff_allele)) %>% 
  # remove "Do/Re" column
  dplyr::select(-`Do/Re`) %>% 
  #eff_allele컬럼에 dbsnp의 ref/alt 확인하여 없으면 strand 변경 
  #### khy
  # for "%dopar% ", need "registerDoParallel(slave.cluster)"
  #registerDoParallel(slave.cluster)
  mutate(eff_allele = foreach(row = 1:n(), .combine = 'c', .multicombine = TRUE) %do% 
                        findAltDbsnp(SNP[row], eff_allele[row])) %>% 
  #mutate(eff_allele = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %dopar% 
  #                      findAltDbsnp(SNP[row], eff_allele[row])) %>% 
  #eff_allele컬럼에 dbsnp에서 없으면 alt sequence로 입력
  #### khy
  # for "%dopar% ", need "registerDoParallel(slave.cluster)"
  #registerDoParallel(slave.cluster)
  mutate(eff_allele = foreach(row = 1:n(), .combine = 'c', .multicombine = TRUE) %do% 
                        checkAlleleDbsnp(SNP[row], eff_allele[row])) %>% 
  #mutate(eff_allele = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %dopar% 
  #                      checkAlleleDbsnp(SNP[row], eff_allele[row])) %>% 
  # if "eff_allele" column is empty, remove that row
  filter(eff_allele != '') %>% 
  # convert to data table
  data.table %>% 
  (function(dat) {
    # set "SNP" key
    setkey(dat, SNP)
    ##### khy ????
    # save for debug
    fwrite(dat, file = "beforebim.csv")
    #fwrite(dat, file = "beforebim.Rdata")
    cat("##### before bim saved", sprintf('%s', Sys.time()), "\n")
    dat
  }) %>%
'''

########################################
print('#### starting to get "unfiltered_meta"\n')

g1 = pd.read_csv(g1_file, encoding="cp949")
print(g1.shape)
print(g1.head())

g1['Bio.Serial'] = np.arange(1, len(g1) + 1)
g1 = g1[g1['code'].notnull()].reset_index(drop=True)
g1 = g1[g1['1st check'] == 'O'].reset_index(drop=True)

# regular expression: '^rs\d+$' vs '^rs\\d+$'
# in R, must use '^rs\\d+$', NOT '^rs\d+$'
# in python, use '^rs\d+$', but '^rs\\d+$' may work
# r grepl => python equivalent
# [ref] https://stackoverflow.com/questions/42853457/python-equivalent-to-grepl-r-with-condition/42853670#42853670
select_snps = g1['RSID'].isnull() & [re.search('^rs\d+\s*$', s) is not None for s in g1['SNPS']]
# select_snps = g1['RSID'].isnull() & [re.search('^rs\d+$', s) is not None for s in g1['SNPS']]
###### above comment out line have bug:
# "SNPS" column may have 'rs1234 ', not 'rs1234',
# in case of 'rs1234 ', "re.search('^rs\d+$', s)" give "FALSE", which is not what we want
# to fix this bug, must consider ending space character,
# so must use "re.search('^rs\d+\s*$', s)"

# for dataframe column selection, g1['RSID'], g1.RSID all work well
g1['RSID'] = np.where(select_snps, g1['SNPS'], g1['RSID'])

rsid_valid = [re.search('rs\d+$', s) is not None for s in g1['RSID']]
g1 = g1[rsid_valid].reset_index(drop=True)

# rename column name, need this in python?
g1 = g1.rename(columns={'RSID' : 'SNP', 'Asia(EAS)' : 'sample.EAS', \
                        'Non-EAS' : 'sample.OTHER', 'Total sample size' : 'sample.TOTAL'})

g1['sample.EAS'] = np.where(g1['sample.EAS'].isnull(), 0, g1['sample.EAS'])
# also try below line
# g1.loc[g1['sample.EAS'].isnull(), 'sample.EAS'] = 0
g1['sample.OTHER'] = np.where(g1['sample.OTHER'].isnull(), 0, g1['sample.OTHER'])
g1['sample.TOTAL'] = np.where(g1['sample.TOTAL'].isnull(), 0, g1['sample.TOTAL'])

##### to check "g1", you can save "g1" into csv file:
# "g1.to_csv('xxx.csv', encoding='cp949', index=False)"

# 'P-VALUE' column may have '' or '2..68E-9'
# "errors='coerce'" coerce bad/missing value to Nan 
# [ref] https://stackoverflow.com/questions/15891038/change-column-type-in-pandas
g1['P-VALUE'] = pd.to_numeric(g1['P-VALUE'], errors='coerce')
g1 = g1[g1['P-VALUE'].notnull()].reset_index(drop=True)

# select column
g1 = g1[['code', 'SNP', 'sample.EAS', 'sample.OTHER', 'sample.TOTAL', \
          'OR/RR/BETA', 'OR or BETA', 'Effect', 'P-VALUE', 'PUBMEDID', \
          'LINK', 'Do/Re', 'Bio.Serial']]

# rename column name: 'Effect' -> 'eff_allele'
g1 = g1.rename(columns={'Effect' : 'eff_allele'})

# pattern string: one of 'beta', '% effect', 'coef', 'AUC', 'delta', 'mean', 'estimate', 'h2', 'u', 'z', 'ln(OR)'
pattern_str = '(beta|% effect|coef|AUC|delta|mean|estimate|h2|u|z|ln\(OR\))'
fill_beta_str = g1['OR/RR/BETA'].notnull() & \
  [re.search(pattern_str, s, re.IGNORECASE) is not None for s in g1['OR/RR/BETA']]
# 're.IGNORECASE' flag in "re.search":
# [ref] https://stackoverflow.com/questions/500864/case-insensitive-regular-expression-without-re-compile

g1['OR/RR/BETA'] = np.where(fill_beta_str, 'BETA', g1['OR/RR/BETA'])

g1['OR or BETA'] = pd.to_numeric(g1['OR or BETA'], errors='coerce')
g1 = g1[g1['OR or BETA'].notnull()].reset_index(drop=True)
print('### before log, g1 length =', len(g1))
g1.to_csv('xx_py.csv', encoding='cp949', index=False)

beta_str_exist = [re.search('beta', s, re.IGNORECASE) is not None for s in g1['OR/RR/BETA']]
# below "notnull" is not necessary because already filtered in above line
# beta_str_exist = g1['OR or BETA'].notnull() & \
#  [re.search('beta', s, re.IGNORECASE) is not None for s in g1['OR/RR/BETA']]
g1['BETA'] = np.where(beta_str_exist, g1['OR or BETA'], np.log(g1['OR or BETA']))
# above line give warning in "np.log": "divide by zero", "invalid value"
# this come when log input is negative: log(-2) => nan

# remove nan in 'BETA' column, which come from negative log input
beta_valid = g1['BETA'].notnull()
g1 = g1[beta_valid].reset_index(drop=True)
print('### after log, g1 length =', len(g1))

# "errors='coerce'" give 'NaN' when 'PUBMEDID' is not digit string
g1['PUBMEDID'] = pd.to_numeric(g1['PUBMEDID'], errors='coerce')
pubmedid_valid = g1['PUBMEDID'].notnull()

for i in range(len(g1)):
  if pubmedid_valid[i]:
    g1.loc[i, 'LINK'] = 'www.ncbi.nlm.nih.gov/pubmed/{:d}'.format(int(g1['PUBMEDID'][i]))

print('### after pubmedid, g1 length =', len(g1))

g1.to_csv('xxx_py.csv', encoding='cp949', index=False)

# "str(s)": need for when "s" is not string
g1['DOMINANT'] = [re.search('Do', str(s), re.IGNORECASE) is not None for s in g1['Do/Re']]

g1['RECESSIVE'] = [re.search('Re', str(s), re.IGNORECASE) is not None for s in g1['Do/Re']]

# sys.exit()

# convert uppercase in 'eff_allele' column
g1['eff_allele'] = g1['eff_allele'].str.upper()

g1.drop(columns= ['Do/Re'])

for i in range(len(g1)):
  g1.loc[i, 'eff_allele'] = findAltDbsnp(g1['SNP'][i], g1['eff_allele'][i])
  # below line give warning:
  # "See the caveats in the documentation: 
  # https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy"
  # g1['eff_allele'][i] = findAltDbsnp(g1['SNP'][i], g1['eff_allele'][i])

for i in range(len(g1)):
  g1.loc[i, 'eff_allele'] = checkAlleleDbsnp(g1['SNP'][i], g1['eff_allele'][i])
  # g1['eff_allele'][i] = checkAlleleDbsnp(g1['SNP'][i], g1['eff_allele'][i])

g1 = g1[g1['eff_allele'].notnull()].reset_index(drop=True)

g1.to_csv('before_bim_py.csv', encoding='utf-8', index=False)
########### [important] difference between two csv file:
# "beforebim.csv" from R code, row length = 856
# "before_bim_py.csv" from python code, row length = 31815
# row data which exist in both csv file(same "code", same "SNP") is same
# why????? i bet on python code

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

#### ref_bim_file not have column name, separation character = "tab", column name prefix = 'V'
bim = pd.read_csv(ref_bim_file, encoding="UTF-8", header=None, sep='\t', prefix='V')
print(bim.shape)
print(bim.head())

# select column from "bim" dataframe: R = one-indexing, python = zero-indexing
bim = bim[['V0', 'V1', 'V3', 'V4', 'V5']]

# rename column name
bim = bim.rename(columns={'V0':'CHR', 'V1':'snpname', 'V3':'POS', 'V4':'REF', 'V5':'ALT'})

