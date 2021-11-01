
import pandas as pd
import numpy as np

from file_db_define import *

'''
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
'''

def mutate_TIMES_column(code_continuous_list, code_continous, BETA, code, LINK):

  return (pd.notna(BETA)) & \
   ( \
      (code not in code_continuous_list) | \
      (sum( \
            (code_continous['code'] == code) & \
            (code_continous['LINK'] == LINK) & \
            ~(pd.isna(code_continous['mean']) | (code_continous['mean'] == 0)) \
          ) > 0 \
      ) \
   )  


'''
  #mutate(BETA = as.numeric(BETA)) %>%
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


def mutate_BETA_column(code_continuous_list, code_continous, bim, i):

  if sum(bim.TIMES[bim.code == bim.code[i]]) == 0:
    return (0.001 * bim.BETA[i] / abs(bim.BETA[i]))

  elif bim.code[i] not in code_continuous_list:
    return bim.BETA[i]

  select_row = (code_continous['code'] == bim.code[i]) & (code_continous['LINK'] == bim.LINK[i])
  mean = code_continous['mean'][select_row].values
  # "code_continous['code']", "code_continous.code]" is all ok,
  # but dont use "code_continous.mean", python recognize "mean" as method
  #### in ipython, if you test and change function which is imported in other file,
  #### must quit ipython, and restart

  if len(mean) == 0:
    return (0.001 * bim.BETA[i] / abs(bim.BETA[i]))

  mean = mean[0]
  if mean == 0:
    return (0.001 * bim.BETA[i] / abs(bim.BETA[i]))
  elif mean < 0:
    return bim.BETA[i]
  else:
    return (np.log((bim.BETA[i] + mean) / mean))


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


####### what "G1_TIMES" is for? "G1_TIMES" is empty. right? 
####### R code produce non empty "G1_TIMES"
def write_G1_TIMES(bim, G1_TIMES_file):

  # in "pd.groupby", meaning of "as_index=False"
  # [ref] https://stackoverflow.com/questions/41236370/what-is-as-index-in-groupby-in-pandas/41237258
  bim = bim.groupby('code', as_index=False)[['TIMES']].count()
  bim = bim[bim['TIMES'] == 0].reset_index(drop=True)
  bim[['code']].to_csv(G1_TIMES_file, encoding='utf-8', index=False)


def mutate_racial_column(bim):

  '''
  mutate(
      East.Asian.Ancestry.Sample.Initial = 
        ifelse(East.Asian.Ancestry.Sample.Initial %in% c(0, NA),
          ifelse(!is.na(sample.EAS), sample.EAS, 0), East.Asian.Ancestry.Sample.Initial),
      East.Asian.Ancestry.Sample.Replication = 
        ifelse(East.Asian.Ancestry.Sample.Replication %in% c(0, NA),
          0, East.Asian.Ancestry.Sample.Replication)) %>%
  '''

  x = bim['East.Asian.Ancestry.Sample.Initial']
  x_cond = (x == 0) | (x.isnull())
  bim['East.Asian.Ancestry.Sample.Initial'] = np.where(x_cond, \
    np.where(bim['sample.EAS'].notnull(), bim['sample.EAS'], 0), \
      bim['East.Asian.Ancestry.Sample.Initial'])

  x = bim['East.Asian.Ancestry.Sample.Replication']
  x_cond = (x == 0) | (x.isnull())
  bim['East.Asian.Ancestry.Sample.Replication'] = np.where(x_cond, \
    0, bim['East.Asian.Ancestry.Sample.Replication'])

  '''
    mutate(
      Asian.Ancestry.Sample.Initial = 
        ifelse(Asian.Ancestry.Sample.Initial %in% c(0, NA),
          ifelse(is.na(East.Asian.Ancestry.Sample.Initial), 0, East.Asian.Ancestry.Sample.Initial),
          Asian.Ancestry.Sample.Initial),
      Asian.Ancestry.Sample.Replication = 
        ifelse(Asian.Ancestry.Sample.Replication %in% c(0, NA),
          ifelse(is.na(East.Asian.Ancestry.Sample.Replication), 
            0, East.Asian.Ancestry.Sample.Replication),
          Asian.Ancestry.Sample.Replication),
      European.Ancestry.Sample.Initial = 
        ifelse(European.Ancestry.Sample.Initial %in% c(0, NA), 
          0, European.Ancestry.Sample.Initial),
      European.Ancestry.Sample.Replication = 
        ifelse(European.Ancestry.Sample.Replication %in% c(0, NA), 
          0, European.Ancestry.Sample.Replication)) %>%
  '''

  x = bim['Asian.Ancestry.Sample.Initial']
  x_cond = (x == 0) | (x.isnull())
  bim['Asian.Ancestry.Sample.Initial'] = np.where(x_cond, \
    np.where(bim['East.Asian.Ancestry.Sample.Initial'].notnull(), bim['East.Asian.Ancestry.Sample.Initial'], 0), \
      bim['Asian.Ancestry.Sample.Initial'])

  x = bim['Asian.Ancestry.Sample.Replication']
  x_cond = (x == 0) | (x.isnull())
  bim['Asian.Ancestry.Sample.Replication'] = np.where(x_cond, \
    np.where(bim['East.Asian.Ancestry.Sample.Replication'].notnull(), bim['East.Asian.Ancestry.Sample.Replication'], 0), \
      bim['Asian.Ancestry.Sample.Replication'])
  
  x = bim['European.Ancestry.Sample.Initial']
  x_cond = (x == 0) | (x.isnull())
  bim['European.Ancestry.Sample.Initial'] = np.where(x_cond, \
    0, bim['European.Ancestry.Sample.Initial'])

  x = bim['European.Ancestry.Sample.Replication']
  x_cond = (x == 0) | (x.isnull())
  bim['European.Ancestry.Sample.Replication'] = np.where(x_cond, \
    0, bim['European.Ancestry.Sample.Replication'])

  '''
    mutate(
      Total.Ancestry.Sample.Initial = 
        ifelse(Total.Ancestry.Sample.Initial %in% c(0, NA),
          ifelse(!is.na(sample.TOTAL), 
            sample.TOTAL, 
            Asian.Ancestry.Sample.Initial + European.Ancestry.Sample.Initial),
          Total.Ancestry.Sample.Initial),
      Total.Ancestry.Sample.Replication = 
        ifelse(Total.Ancestry.Sample.Replication %in% c(0, NA),
          Asian.Ancestry.Sample.Replication + European.Ancestry.Sample.Replication,
          Total.Ancestry.Sample.Replication)) %>%
  '''

  x = bim['Total.Ancestry.Sample.Initial']
  x_cond = (x == 0) | (x.isnull())
  bim['Total.Ancestry.Sample.Initial'] = np.where(x_cond, \
    np.where(bim['sample.TOTAL'].notnull(), bim['sample.TOTAL'], bim['Asian.Ancestry.Sample.Initial'] + bim['European.Ancestry.Sample.Initial']), \
      bim['Total.Ancestry.Sample.Initial'])

  x = bim['Total.Ancestry.Sample.Replication']
  x_cond = (x == 0) | (x.isnull())
  bim['Total.Ancestry.Sample.Replication'] = np.where(x_cond, \
    bim['Asian.Ancestry.Sample.Replication'] + bim['European.Ancestry.Sample.Replication'], bim['Total.Ancestry.Sample.Replication'])

  '''
    #### khy
    # study: (East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication > 0) * 15
    mutate(Point = (
            ifelse(is.na(East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication), 
            0, (East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication > 0) * 15) +
            ifelse(is.na(Asian.Ancestry.Sample.Initial + Asian.Ancestry.Sample.Replication), 
            0, (Asian.Ancestry.Sample.Initial + Asian.Ancestry.Sample.Replication > 0) * 30) +
            ifelse(is.na(European.Ancestry.Sample.Initial + European.Ancestry.Sample.Replication), 
            0, (European.Ancestry.Sample.Initial + European.Ancestry.Sample.Replication > 0) * 10) +
            ifelse(is.na(Total.Ancestry.Sample.Initial + Total.Ancestry.Sample.Replication), 
            0, pmax(pmin(log10(Total.Ancestry.Sample.Initial + Total.Ancestry.Sample.Replication), 5), 0) * 10) +
            ifelse(is.na(`P-VALUE`), 
            0, pmin(-log10(`P-VALUE`), 10))) / 1.15) %>%
  '''

  cond_1 = (bim['East.Asian.Ancestry.Sample.Initial'] + bim['East.Asian.Ancestry.Sample.Replication']).isnull()
  # python operator precedence is same as R: (high) +, (low) >
  point_1 = np.where(cond_1, \
    0, (bim['East.Asian.Ancestry.Sample.Initial'] + bim['East.Asian.Ancestry.Sample.Replication'] > 0) * 15)

  cond_2 = (bim['Asian.Ancestry.Sample.Initial'] + bim['Asian.Ancestry.Sample.Replication']).isnull()
  point_2 = np.where(cond_2, \
    0, (bim['Asian.Ancestry.Sample.Initial'] + bim['Asian.Ancestry.Sample.Replication'] > 0) * 30)

  cond_3 = (bim['European.Ancestry.Sample.Initial'] + bim['European.Ancestry.Sample.Replication']).isnull()
  point_3 = np.where(cond_3, \
    0, (bim['European.Ancestry.Sample.Initial'] + bim['European.Ancestry.Sample.Replication'] > 0) * 10)

  pmin = np.minimum(np.log10(bim['Total.Ancestry.Sample.Initial'] + bim['Total.Ancestry.Sample.Replication']), 5)
  cond_4 = (bim['Total.Ancestry.Sample.Initial'] + bim['Total.Ancestry.Sample.Replication']).isnull()
  point_4 = np.where(cond_4, 0, np.maximum(pmin, 0) * 10)

  point_5 = np.where(bim['P-VALUE'].isnull(), 0, np.minimum(-np.log10(bim['P-VALUE']), 10))

  bim['Point'] = (point_1 + point_2 + point_3 + point_4 + point_5) / 1.15

  return bim

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

def mutate_beta_zero_column(bim, dat_ld):

  unfiltered_meta = pd.DataFrame()
  CODE = bim['code'].unique()
  CODE_len = len(CODE)
  for count, code in enumerate(CODE):
    print('#### code = {}, {} / {}'.format(code, count, CODE_len))
    bim_code = bim[bim['code'] == code]
    bim_code = bim_code.sort_values(by=['Point', 'P-VALUE'], ascending=[False, False]) \
      .reset_index(drop=True)
    ###############################################################################################
    ########### MUST USE above ".reset_index(drop=True)"
    ########### only when ".reset_index(drop=True)" is used, "loc" may give what you want
    ########### i'm beginner in pandas: "loc", "iloc" is very different!!!!
    ########### for non-monotonic integer index, "loc" give key error
    ########### [ref] https://stackoverflow.com/questions/31593201/how-are-iloc-and-loc-different
    ########### MUST REVIEW other code where "loc" was used:
    ########### "process_bim.py", "process_g1.py", "get_code_continuous_mean.py", ...
    #####################################################################################
    for i in range(len(bim_code)):
      # print(i)
      if i == 0:
        bim_code.loc[i, 'beta_zero'] = False
      else:
        snpname_list = list(bim_code['snpname'][0:i])
        # print(snpname_list)
        snp_a_match = dat_ld['SNP_A'].isin(snpname_list)
        snp_b_list = list(dat_ld[snp_a_match]['SNP_B'])
        # print(snp_b_list)
        # print(bim_code['snpname'][i])
        ######### very long run time, last two line with zero was printed (why?)
        bim_code.loc[i, 'beta_zero'] = (bim_code['snpname'][i] in snp_b_list)

    bim_code['BETA'] = np.where(bim_code['beta_zero'], 0, bim_code['BETA'])
    bim_code['SNP_OR'] = np.where(bim_code['beta_zero'], 1, bim_code['SNP_OR'])

    # [ref] https://pandas.pydata.org/docs/reference/api/pandas.concat.html
    # R: rbindlist(list(...)), python: df3 = pd.concat([df1, df2], ignore_index=True)
    unfiltered_meta = pd.concat([unfiltered_meta, bim_code], ignore_index=True)

  return unfiltered_meta

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

def filter_mutate(group):

  '''
  filter(!duplicated(snpname)) %>%
  mutate(order_id = 1:n()) %>%
  '''

  group = group.drop_duplicates(subset=['snpname'])
  group['order_id'] = np.arange(1, len(group) + 1)

  return group


'''
jointPAF <- function(PAFs) {
  rowSums(PAFs / (1 - PAFs)) / (1 + rowSums(PAFs / (1 - PAFs)))
}
'''

def jointPAF(PAFs):

  paf = PAFs / (1 - PAFs)
  snp_paf = paf.sum(axis=1) / (1 + paf.sum(axis=1))

  return snp_paf


'''
  mutate(snp_paf = 
    jointPAF(cbind(ifelse(ab_or > 1, 
                    2 * snp_freq * (1 - snp_freq) * (ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (ab_or - 1)), 
                    2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1))),
                  ifelse(bb_or > 1, 
                    snp_freq * snp_freq * (bb_or - 1) / (1 + snp_freq * snp_freq * (bb_or - 1)), 
                    (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1) / (1 + (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1)))
                  )
            )
        )
'''


def mutate_snp_paf(unfiltered_meta):

  snp_freq = unfiltered_meta['snp_freq'].values
  ab_or = unfiltered_meta['ab_or'].values
  bb_or = unfiltered_meta['bb_or'].values

  x = 2 * snp_freq * (1 - snp_freq) * (ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (ab_or - 1))
  y = 2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1))
  column_ab = pd.DataFrame(np.where(ab_or > 1, x, y))

  x = snp_freq * snp_freq * (bb_or - 1) / (1 + snp_freq * snp_freq * (bb_or - 1))
  y = (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1) / (1 + (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1))
  column_bb = pd.DataFrame(np.where(bb_or > 1, x, y))

  column = pd.concat([column_ab, column_bb], axis=1)

  # 'snp_paf' column have 2 nan
  unfiltered_meta['snp_paf'] = jointPAF(column)

  return unfiltered_meta


if __name__ == "__main__":

  test_mutate_TIMES_column = False
  test_unfiltered_meta_group_by = False
  test_mutate_snp_paf = True

  if test_mutate_TIMES_column:

    code_continous = pd.read_csv(code_continous_file, encoding="UTF-8")
    print(code_continous.shape)
    print(code_continous.head())

    code_continuous_list = pd.read_csv(continuous_file, encoding="UTF-8").code.unique()
    # "unique()" return numpy array
    # code_continuous_list = pd.DataFrame(code_continuous_list)
    print(code_continuous_list.shape)
    # print(code_continuous_list.head())

    BETA = 0.01
    code = 'CGH'
    LINK = 'www.ncbi.nlm.nih.gov/pubmed/29403010'
    TIMES = mutate_TIMES_column(code_continuous_list, code_continous, BETA, code, LINK)
    print(TIMES)

  if test_unfiltered_meta_group_by:

    unfiltered_meta = pd.read_csv(after_mutate_beta_zero_file, encoding='utf-8')

    '''
    group_by(code) %>%
    filter(!duplicated(snpname)) %>%
    mutate(order_id = 1:n()) %>%
    ungroup %>%
    mutate(gene_description = 'Not available') %>%
    mutate(title = 'NA') %>%
    rename(description = gene_description) %>%
    '''

    unfiltered_meta = unfiltered_meta.groupby('code', as_index=False).apply(filter_mutate).reset_index(drop=True)

  if test_mutate_snp_paf:

    unfiltered_meta = pd.read_csv(before_mutate_snp_paf_file, encoding='utf-8')

    unfiltered_meta = mutate_snp_paf(unfiltered_meta)




