
import sys
import pandas as pd

from file_db_define import *

def get_code_continuous_mean(g1_continuous_mean_file):

  '''
  # process code continuous mean value
  code.continuous <- fread(G1_mean_csv_file, encoding = "UTF-8") %>%
    (function(dat) {
      dat[, .(code, PUBMEDID, LINK, mean, `DISEASE/TRAIT`)] %>%
        merge(dat %>%
                ############# khy ???
                # comment out because it make program malfunction
                # why? suspect "grepl" wrong
                # but in R terminal, it seems work well
                filter(!grepl('부적합', `Comment`)) %>%
                group_by(code) %>%
                summarize(n = sum(!is.na(mean))) %>%
                filter(n > 0), 
                by = c('code'), all.x = FALSE, all.y = TRUE)
    }) %>%
    filter(!is.na(mean)) %>%
    mutate(LINK = ifelse(is.na(LINK), '', LINK)) %>%
    filter(!is.na(code)) %>%
    mutate(LINK = ifelse(!is.na(PUBMEDID) & PUBMEDID != '' & PUBMEDID != '0', 
      sprintf('www.ncbi.nlm.nih.gov/pubmed/%s', PUBMEDID), LINK))

  save(code.continuous, file = 'code.continuous.RData')
  '''

  g1_continuous_mean = pd.read_csv(g1_continuous_mean_file, encoding="cp949")
  print(g1_continuous_mean.shape)
  print(g1_continuous_mean.head())

  x = g1_continuous_mean[['code', 'PUBMEDID', 'LINK', 'mean', 'DISEASE/TRAIT']]
  print(x.shape)
  print(x.head())

  y = g1_continuous_mean[g1_continuous_mean['Comment'] != '결과 부적합']
  print(y.shape)
  # in "pd.groupby", meaning of "as_index=False"
  # [ref] https://stackoverflow.com/questions/41236370/what-is-as-index-in-groupby-in-pandas/41237258
  y = y.groupby('code', as_index=False)[['mean']].count()

  # need "reset_index(drop=True)"? absolutely yes
  # rename column name: 'mean' -> 'n'
  y = y[y['mean'] > 0].rename(columns={'mean' : 'n'}).reset_index(drop=True)
  # y = y.loc[y['mean'] > 0].rename(columns={'mean' : 'n'}).reset_index(drop=True)

  x = x.merge(y, on='code', how='right')

  x = x[x['mean'].notnull()].reset_index(drop=True)
  # '~' tilde operator in pandas, [ref] https://blog.finxter.com/tilde-python/
  # x = x[~x['mean'].isnull()]

  # [ref] https://towardsdatascience.com/python-pandas-vs-r-dplyr-5b5081945ccb
  x.loc[x['LINK'].isnull(), 'LINK'] = ''

  x = x[x['code'].notnull()].reset_index(drop=True)
  # x = x[~x['code'].isnull()]

  # sys.exit()

  '''
  # "errors='coerce'" give 'NaN' when 'PUBMEDID' is not digit string
  g1['PUBMEDID'] = pd.to_numeric(g1['PUBMEDID'], errors='coerce')
  pubmedid_valid = g1['PUBMEDID'].notnull()
  g1 = g1[pubmedid_valid].reset_index(drop=True)
  '''

  p_valid = x['PUBMEDID'].notnull() & (x['PUBMEDID'] != 0)
  x = x[p_valid].reset_index(drop=True)

  '''
  # 'p': pandas series, 'p' data type is 'float64', not stirng
  p = x['PUBMEDID'].reset_index(drop=True)

  p_valid = (p.notnull() & (p != 0)).reset_index(drop=True)
  '''

  # i know what to do, but dont know how to do
  # [solved] i knew how to do
  # after merge: "x = x.merge(y, on='code', how='right')", 
  # x['PUBMEDID'] data type was changed from string to float, so need to cast type to "int"
  for i in range(len(x)):
      x.loc[i, 'LINK'] = 'www.ncbi.nlm.nih.gov/pubmed/{:d}'.format(int(x['PUBMEDID'][i]))

  '''
  for i in range(len(x)):
    if p_valid[i]:
      x.loc[i, 'LINK'] = 'www.ncbi.nlm.nih.gov/pubmed/{:d}'.format(int(p[i]))
  '''

  # x.to_csv(code_continous_file, index=False)

  return x


if __name__ == "__main__":

  # from file_db_define import *

  # g1_continuous_mean_file = "C:/Users/hykang/Desktop/R_study/reference/fixed_G1_Continuous_Mean.csv"
  # code_continous_file = "C:/Users/hykang/Desktop/R_study/code_continuous.csv"

  code_continuous = get_code_continuous_mean(g1_continuous_mean_file)

  code_continuous.to_csv(code_continous_file, index=False)
