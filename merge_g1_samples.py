
import pandas as pd
import numpy as np

from check_LD_Plink import *
from mutate_foreach import *
from file_db_define import *
from plink_util import *

def merge_g1_samples(bim, g1_samples_file):

  '''
    merge(fread('reference/g1_samples.csv', encoding = "UTF-8") %>% dplyr::select(-race, -Outcomes),
          by = c('PUBMEDID', 'LINK', 'code'), all.x = TRUE, all.y = FALSE)
    }) %>%
  '''

  g1_smpl = pd.read_csv(g1_samples_file, encoding="UTF-8")
  g1_smpl = g1_smpl.drop(columns=['race', 'Outcomes'])
  print(g1_smpl.shape)
  print(g1_smpl.head())

  # merge "bim" with "g1_smpl": on = 'SNP', how = 'right'('right' is right?)
  # how option = ‘left’, ‘right’, ‘outer’, ‘inner’
  # [ref] https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html
  # [ref] https://stackoverflow.com/questions/53645882/pandas-merging-101/53645883#53645883
  bim = bim.merge(g1_smpl, on=['PUBMEDID', 'LINK', 'code'], how='left')

  bim.to_csv(after_merge_g1_samples_file, encoding='utf-8', index=False)

  bim = mutate_racial_column(bim)

  '''
    (function(dat) {
    ####### khy
    # debug
    fwrite(dat, file = "before_checkLDPlink.csv")
    #fwrite(dat, file = "before_checkLDPlink.rdata")
    cat("######## before checkLDPlink call saved", sprintf('%s', Sys.time()), "\n")

    dat.ld <- checkLDPlink(dat$snpname %>% unique %>% setdiff(c(NA, '')))
    ##### khy
    # debug
    fwrite(dat.ld, file = "checkLDPlink_output.csv")
    # wait for fwrite
    #Sys.sleep(3)
    cat("######## checkLDPlink output saved", sprintf('%s', Sys.time()), "\n")
  '''

  bim.to_csv(before_checkLDPlink_file, encoding='utf-8', index=False)

  # extract "snpname" column, and make unique
  snpname = bim['snpname'].unique()

  # above "test_input" returned data type is "numpy.ndarray", 
  # so convert to dataframe, and remove NA
  snp = pd.DataFrame(snpname).dropna()

  dat_ld = check_LD_Plink(snp)

  dat_ld.to_csv(checkLDPlink_output_file, encoding='utf-8', index=False)

  '''
      foreach(CODE = unique(dat$code), 
            .combine = function(...) rbindlist(list(...)), .multicombine=TRUE) %do% {
      dat.code <- dat %>%
        filter(code == CODE)
      dat.code %>%
        arrange(desc(Point), desc(`P-VALUE`)) %>%
        mutate(beta_zero = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %do% {
          if (row == 1)
            return(FALSE)
          return(snpname[row] %in% dat.ld[SNP_A %in% snpname[1:(row - 1)], SNP_B])
        }) %>%
        mutate(BETA = ifelse(beta_zero, 0, BETA)) %>%
        mutate(SNP_OR = ifelse(beta_zero, 1, SNP_OR))
    }
  }) %>%
  '''

  unfiltered_meta = mutate_beta_zero_column(bim, dat_ld)
  # "mutate_beta_zero_column" run long time
  # so save "unfiltered_meta" not to wait for "mutate_beta_zero_column" to finish
  unfiltered_meta.to_csv(after_mutate_beta_zero_file, encoding='utf-8', index=False)

  '''
  group_by(code) %>%
  filter(!duplicated(snpname)) %>%
  mutate(order_id = 1:n()) %>%
  ungroup %>%
  mutate(gene_description = 'Not available') %>%
  mutate(title = 'NA') %>%
  rename(description = gene_description) %>%
  '''

  # in "pd.groupby('code', as_index=False)", meaning of "as_index=False":
  # when "as_index=False", index is monotonic integer,
  # when "as_index=True", index is 'code'.
  # [ref] https://stackoverflow.com/questions/41236370/what-is-as-index-in-groupby-in-pandas/41237258
  unfiltered_meta = unfiltered_meta.groupby('code', as_index=False) \
    .apply(filter_mutate).reset_index(drop=True)
  # ".apply(filter_mutate)": [ref] https://johnpaton.net/posts/groupby-without-aggregation/
  # ".reset_index(drop=True)" is same as "ungroup" in R
  # [ref] https://stackoverflow.com/questions/20122521/is-there-an-ungroup-by-operation-opposite-to-groupby-in-pandas

  unfiltered_meta.to_csv(after_groupby_unfiltered_meta_file, encoding='utf-8', index=False)

  unfiltered_meta['gene_description'] = 'Not available'
  unfiltered_meta['title'] = 'NA'
  unfiltered_meta = unfiltered_meta.rename(columns={'gene_description' : 'description'})

  '''
  mutate(snp_freq = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %do% {
    con <- dbConnect(SQLite(), 'reference/Plink/ref.frq.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    res <- dbGetQuery(con, 
              sprintf('select * from freq where SNP="%s" and ((A1="%s" and A2="%s") or (A1="%s" and A2="%s"))',
              snpname[row], ref_allele[row], eff_allele[row], eff_allele[row], ref_allele[row]))
    if (nrow(res) == 0)
      return(NA)
    else {
      res <- res[1,]
      return(ifelse(res$A1 == eff_allele[row], res$MAF, 1 - res$MAF))
    }
  }) %>%
  '''

  for i in range(len(unfiltered_meta)):
    snpname = unfiltered_meta['snpname'][i]
    ref_allele = unfiltered_meta['ref_allele'][i]
    eff_allele = unfiltered_meta['eff_allele'][i]
    unfiltered_meta.loc[i, 'snp_freq'] = \
      query_plink_ref_freq(plink_ref_freq_db, snpname, ref_allele, eff_allele)

  '''
  rename(snp_or = SNP_OR) %>%
  mutate(snp_freq = pmin(snp_freq, 0.99999)) %>%
  mutate(snp_freq = pmax(snp_freq, 0.00001)) %>%
  mutate(aa_or = 1) %>%
  mutate(ab_or = ifelse(RECESSIVE, 1, snp_or)) %>%
  mutate(bb_or = ifelse(DOMINANT | RECESSIVE, snp_or, snp_or**2)) %>%
  mutate(lang = 'ko') %>%
  mutate(version = 1) %>%
  rename(ref_url = LINK) %>%
  '''

  unfiltered_meta = unfiltered_meta.rename(columns={'SNP_OR' : 'snp_or'})
  unfiltered_meta['snp_freq'] = np.minimum(unfiltered_meta['snp_freq'], 0.99999)
  unfiltered_meta['snp_freq'] = np.maximum(unfiltered_meta['snp_freq'], 0.00001)
  unfiltered_meta['aa_or'] = 1
  unfiltered_meta['ab_or'] = np.where(unfiltered_meta['RECESSIVE'], 1, unfiltered_meta['snp_or'])
  unfiltered_meta['bb_or'] = np.where(unfiltered_meta['RECESSIVE'] | unfiltered_meta['RECESSIVE'], \
    unfiltered_meta['snp_or'], unfiltered_meta['snp_or']**2)
  unfiltered_meta['lang'] = 'ko'
  unfiltered_meta['version'] = 1
  unfiltered_meta = unfiltered_meta.rename(columns={'LINK' : 'ref_url'})

  unfiltered_meta.to_csv(before_mutate_snp_paf_file, encoding='utf-8', index=False)

  unfiltered_meta = mutate_snp_paf(unfiltered_meta)

  unfiltered_meta.to_csv(unfiltered_meta_file, encoding='utf-8', index=False)

  return unfiltered_meta


if __name__ == "__main__":

  bim = pd.read_csv(before_merge_g1_samples_file, encoding="UTF-8")
  print(bim.shape)
  print(bim.head())

  unfiltered_meta = merge_g1_samples(bim, g1_samples_file)


