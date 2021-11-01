
import pandas as pd
import numpy as np
from file_db_define import *

'''
beta_zero <- subset(unfiltered_meta, unfiltered_meta$beta_zero == TRUE)
#### khy, debug
fwrite(beta_zero, file = "beta_zero.csv")
cat("##### beta zero saved", sprintf("%s", Sys.time()), "\n")

final_meta <- subset(unfiltered_meta, unfiltered_meta$beta_zero == FALSE) %>%
  (function(dat) {
    #### khy
    # "%dopar%" -> "%do%"
    foreach(CODE = unique(dat$code), .combine = function(...) rbindlist(list(...)), 
      .multicombine = F) %do% {
      for (cutoff in seq(0.3, 1, 0.02)) {
        res <- dat %>%
          filter(code == CODE) %>%
          filter(snp_paf <= cutoff)
        if (nrow(res) > 1)
          break
      }
      return(res)
    }
  }) %>% bind_rows(beta_zero)
'''


def get_final_meta(unfiltered_meta):

   beta_zero = unfiltered_meta[unfiltered_meta['beta_zero'] == True].reset_index(drop=True)
   beta_zero.to_csv(beta_zero_file, encoding='utf-8', index=False)

   dat = unfiltered_meta[unfiltered_meta['beta_zero'] == False].reset_index(drop=True)

   tmp_meta = pd.DataFrame()
   CODE = dat['code'].unique()
   CODE_len = len(CODE)
   for count, code in enumerate(CODE):
      print('#### code = {}, {} / {}'.format(code, count, CODE_len))
      for cutoff in np.arange(0.3, 1 + 0.02, 0.02):
         dat_filter = dat[dat['code'] == code]
         res = dat_filter[dat_filter['snp_paf'] <= cutoff]
         if len(res) > 1:
            break
  
      tmp_meta = pd.concat([tmp_meta, res], ignore_index=True)

   final_meta = pd.concat([tmp_meta, beta_zero]).reset_index(drop=True)

   final_meta.to_csv(final_meta_file, encoding='utf-8', index=False)

   return final_meta, beta_zero


if __name__ == "__main__":   

   unfiltered_meta = pd.read_csv(unfiltered_meta_file, encoding='utf-8')

   final_meta, beta_zero = get_final_meta(unfiltered_meta)
