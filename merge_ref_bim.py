
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

import pandas as pd
import numpy as np
# re: regular expression
import re, sys

from cytogenic_band import *
from file_db_define import *

def merge_ref_bim(ref_bim_file, g1):

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

  bim.to_csv(afterbim_file, encoding='utf-8', index=False)

  # print(bim.head())

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

  return bim


def mutate_from_cytogenic_band(bim, cytogenic_band_length_file):

  # import pandas as pd
  '''
  from cytogenic_band import \
    location_from_cytogenic_band, ARM_from_cytogenic_band, \
      position_percentage_from_cytogenic_band
  '''

  cytogenic_band_length = pd.read_csv(cytogenic_band_length_file, encoding="UTF-8")
  print(cytogenic_band_length.shape)
  print(cytogenic_band_length.head())

  '''
   (function(dat) {
      #registerDoParallel(detectCores()) 
      ##### khy
      fwrite(dat, file = "ungroup_afterbim.csv")
      #fwrite(dat, file = "ungroup_afterbim.rdata")
      cat("##### ungroup afterbim saved", sprintf('%s', Sys.time()), "\n")
      ### khy
      # %dopar% -> %do%
      ##### khy
      # for ".N" in data.table,
      # see "https://www.rdocumentation.org/packages/data.table/versions/1.10.0/topics/special-symbols"
      dat[, location := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %do% { 
         # if "CHR" not exist in cytogenic band DB, return with "-"
         if (!CHR[row] %in% unique(cytogenic.band.length$chr))
         return('-')
         # else ("CHR" exist in cytogenic band DB)
         row.band <- cytogenic.band.length %>% 
         # select row from "cytogenic.band.length" whose "chr" column is same as "CHR" column
         filter(chr == CHR[row]) %>%  
         # select row whose "bp.start" column is less than "POS" column
         filter(bp.start <= POS[row]) %>% 
         # select row whose "bp.end" column is greater than "POS" column, and save result into row.band
         filter(bp.stop >= POS[row])
         # if row.band is empty,
         if (nrow(row.band) == 0)
         # return with "-"
         return('-')
         # else return with "CHR""arm""band" (example = "4q31.23")
         return(sprintf('%s%s%s', CHR[row], row.band$arm[1], row.band$band[1]))
      }]
      ##### khy
      #registerDoParallel(slave.cluster)
      # %dopar% -> %do%
      dat[, ARM := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %do% { 
         if (!CHR[row] %in% unique(cytogenic.band.length$chr))
         return('-')
         row.band <- cytogenic.band.length %>%
         filter(chr == CHR[row]) %>%
         filter(bp.start <= POS[row]) %>%
         filter(bp.stop >= POS[row])
         if (nrow(row.band) == 0)
         return('-')
         return(row.band$arm[1])
      }]
      #registerDoParallel(slave.cluster)
      ###### khy
      # %dopar% -> %do%
      dat[, position_percentage := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %do% {
         if (!CHR[row] %in% unique(cytogenic.band.length$chr))
         return(0)
         row.band <- cytogenic.band.length %>%
         filter(chr == CHR[row]) %>%
         filter(arm == ARM[row])
         if (nrow(row.band) == 0)
         return(0)
         return((POS[row] - min(row.band$bp.start)) /
                  (max(row.band$bp.stop) - min(row.band$bp.start)))
      }]
      ##### khy
      fwrite(dat, file = "cytogenic_refered.csv")
      #fwrite(dat, file = "cytogenic_refered.rdata")
      cat("###### cytogenic refered saved", sprintf('%s', Sys.time()), "\n")
      dat
      #quit(save = "no")
   })  %>%
  '''

  for i in range(len(bim)):

    bim.loc[i, 'location'] = \
      location_from_cytogenic_band(cytogenic_band_length, bim['CHR'][i], bim['POS'][i])

    bim.loc[i, 'ARM'] = \
      ARM_from_cytogenic_band(cytogenic_band_length, bim['CHR'][i], bim['POS'][i])

    bim.loc[i, 'position_percentage'] = \
      position_percentage_from_cytogenic_band(cytogenic_band_length, \
         bim['CHR'][i], bim['ARM'][i], bim['POS'][i])

  bim.to_csv(cytogenic_refered_file, encoding='utf-8', index=False)

  return bim


if __name__ == "__main__":

  # read "g1" from temporary file
  g1 = pd.read_csv(before_bim_file, encoding="UTF-8")
  print(g1.shape)
  print(g1.head())

  bim = merge_ref_bim(ref_bim_file, g1)

  print("### mutate_from_cytogenic_band")
  bim = mutate_from_cytogenic_band(bim, cytogenic_band_length_file)

