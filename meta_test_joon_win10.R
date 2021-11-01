##### khy
# to run in windows 10
#
# required input file or data:
# (1) 'cytogenic.band.length.RData'
# (2) '/Users/genoplan/metadata/reference/dbsnp/dbsnp.grch37.sqlite'
# (3) 'G1_Continuous_Mean.csv'
# (4) '/Users/genoplan/metadata/reference/Plink/ref.frq'
# (5) '/Users/genoplan/metadata/reference/Plink/ref.frq.sqlite'
# (6) 'code.continuous.RData'
# (7) 'continuous.csv'
# (8) 'g1.csv'
# (9) '/Users/genoplan/metadata/reference/Plink/ref.bim'
# (10) 'g1_samples.csv'
#
# required toolset:
# (1) plink(version 1.9)

#### khy
# for windows os, use not "backslash", but "slash" in "setwd"
#setwd("C:\Users\hykang\Desktop\R study")
setwd("C:/Users/hykang/Desktop/R_study")
cat("## current working directory =", getwd(), "\n")

# /data/98.reference/11.meta_ref/Plink -> /Users/genoplan/metadata/reference/Plink
# /data/98.reference/09.dbsnp/ -> /Users/genoplan/metadata/reference/dbsnp/
# plink -> /Users/genoplan/tools/plink_mac_20210606/plink --bfile

#### khy
# for R markdown: install.packages("rmarkdown")
knitr::opts_chunk$set(echo = TRUE)

#### khy
# must comment out. 
# if not, "readLines('reference/G1_Continuous_Mean.csv')" give error message.
#options(encoding = "UTF-8")

# locale is already set to "Korean". to see, use "Sys.getlocale"
# for R 4.x, 
# use "Sys.setlocale("LC_ALL", "Korean")", "Sys.setlocale("LC_CTYPE", "Korean")"
#Sys.setlocale("LC_ALL", "ko_KR.UTF-8")
#Sys.setlocale("LC_CTYPE", "ko_KR.UTF-8")
# Sys.setenv(LANG = "ko_KR.UTF-8")

#### khy
# for "fread", "fwrite"
#if (!require(data.table)) install.packages("data.table", repos = "https://Rdatatable.gitlab.io/data.table")
#### khy
# delete "repos", why "repos" is needed?
if (!require(data.table)) install.packages("data.table")
library(data.table)

#### khy
# for "mutate", "%>%(pipeline)"
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

### khy
# for "dbConnect", "dbDisconnect", "dbGetQuery", "dbWriteTable"
if (!require(RSQLite)) install.packages('RSQLite')
library(RSQLite)

### khy
# for "registerDoParallel"
if (!require(doParallel)) install.packages('doParallel')
library(doParallel)

if (!require(openxlsx)) install.packages('openxlsx')
library(openxlsx)

if (!require(readr)) install.packages('readr')
library(readr)

##### khy
# for "check.numeric"
if (!require(varhandle)) install.packages('varhandle')
library(varhandle)

##### khy
# for "printf"
#if (!require(R.utils)) install.packages('R.utils')
#library(R.utils)

start_time <- Sys.time()
cat(sprintf("###### program start = %s\n", start_time))

load('reference/cytogenic.band.length.RData')

##### khy
# to use "registerDoParallel":
# (1) install.packages("doParallel")
# (2) library(doParallel)

#registerDoParallel(detectCores())
#registerDoParallel(15)
#### khy
# [ref] https://hoontaeklee.github.io/en/posts/20200607_r%EB%B3%91%EB%A0%AC%EC%B2%98%EB%A6%AC%EB%A9%94%EB%AA%A8/
# NOT USE parallel processing
#n.cores = detectCores() - 1
#slave.cluster = makeCluster(n.cores)
#registerDoParallel(slave.cluster)

##### khy
# to use "foreach":
# (1) install.packages("foreach")
# (2) library(foreach)

# rev.strand
rev.strand <- function(x) {
  foreach(x1 = strsplit(x, '')[[1]], .combine = c) %do% {
    list(A = 'T', T = 'A', C = 'G', G = 'C', D = 'D', I = 'I', '-' = '-')[[x1]]
  } %>%
    paste(collapse = '')
}

# jointPAF
jointPAF <- function(PAFs) {
  rowSums(PAFs / (1 - PAFs)) / (1 + rowSums(PAFs / (1 - PAFs)))
}

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

# checkAlleleDbsnp
checkAlleleDbsnp <- function(SNP, Allele) {
  con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.grch37.sqlite')
  on.exit(dbDisconnect(con), add = TRUE)
  res <- dbGetQuery(con, sprintf('select REF, ALT from dbsnp where `NAME`="%s"', SNP)) %>%
    unlist
  if (Allele %in% res)
    return(Allele)
  else if (rev.strand(Allele) %in% res)
    return(rev.strand(Allele))
  else
    return('')
}

#### khy
# to use "dbConnect", "dbDisconnect", "dbGetQuery":
# (1) install.packages('RSQLite')
# (2) library(RSQLite)

# findAltDbsnp
findAltDbsnp <- function(SNP, Allele) {
  if (grepl('^[ACGT]$', Allele, ignore.case = TRUE))
    return(Allele)
  con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.grch37.sqlite')
  on.exit(dbDisconnect(con), add = TRUE)
  res <- dbGetQuery(con, sprintf('select ALT from dbsnp where `NAME`="%s" limit 1', SNP)) %>%
    unlist
  if (length(res) > 0)
    return(res)
  else
    return(Allele)
}

#### khy
# to use "mutate", "%>%(pipeline)":
# (1) install.packages('dplyr')
# (2) library(dplyr)

#### khy
# "G1_Continuous_Mean.csv" file problem:
# (1) many rows have "Korean", "PUBMEDID", "LINK" data squeezed into "Korean" column
# (2) data in "Korean" column is not displayed correctly

#G1_mean_csv_file <- 'reference/G1_Continuous_Mean.csv'
# use fixed G1 mean csv file:
# to create it, run "fix_G1_mean_csv_file.R" script
G1_mean_csv_file <- 'reference/fixed_G1_Continuous_Mean.csv'

# process code continuous mean value
code.continuous <- fread(G1_mean_csv_file, encoding = "UTF-8") %>%
  (function(dat) {
    dat[, .(code, PUBMEDID, LINK, mean, `DISEASE/TRAIT`)] %>%
      merge(dat %>%
              ############# khy ???
              # comment out because it make program malfunction
              # why? suspect "grepl" wrong
              # but in R terminal, it seems work well
              #filter(!grepl('부적합', `Comment`)) %>%
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
fwrite(code.continuous, file = 'code.continuous.csv')
cat("###### \"code.continuous\" saved", sprintf('%s', Sys.time()), "\n")

# Reference frequencey를 sqlite로 변환

##### khy
# if nay, remove old db
unlink('reference/Plink/ref.frq.sqlite')
fread('reference/Plink/ref.frq', encoding = "UTF-8") %>%
  (function(dat) {
    #### khy
    # db is created
    con <- dbConnect(SQLite(), 'reference/Plink/ref.frq.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    dat %>%
      dbWriteTable(conn = con, name = 'freq', overwrite = TRUE)
    # execute SQL statement: 
    # "CREATE INDEX index_name ON table_name(column_list)"
    dbExecute(conn = con, 'create index idx_freq_SNP_A1_A2 on freq (SNP, A1, A2)')
  })

cat("###### end of sql", sprintf('%s', Sys.time()), "\n")
#quit(save = "no")

# process bio team data
#### khy
# 'code.continuous.RData' was created
# see "save(code.continuous, file = 'code.continuous.RData')"
load('code.continuous.RData')
code.continuous.list <- fread('reference/continuous.csv', encoding = "UTF-8")$code %>% unique

##### khy
# 'reference/g1.csv' file have garbage data, must check

unfiltered_meta <- fread('reference/fixed_g1.csv') %>%
  # Bio.Serial column: Assign serial Number for each SNP
  mutate(Bio.Serial = 1:n()) %>%
  filter(!is.na(code)) %>%
  # Checking history from BS : 1st check is '0'
  filter(`1st check` == 'O') %>%
  
  #### khy
  # what is "RSID %in% c('')"?
  # Substitute RSID col. with SNPS col., if RSID col. is empty and SNPS col. starts with 'rs' and ends with some digits    
  mutate(RSID = ifelse(RSID %in% c('') & grepl('^rs\\d+\\s*$', SNPS), SNPS, RSID)) %>%
  # mutate(RSID = ifelse(RSID %in% c('') & grepl('^rs\\d+$', SNPS), SNPS, RSID)) %>%
  ###### above comment out line have bug:
  # "SNPS" column may have 'rs1234 ', not 'rs1234',
  # in case of 'rs1234 ', "grepl('^rs\\d+$', SNPS)" give "FALSE", which is not what we want
  # to fix this bug, must consider ending space character,
  # so must use "grepl('^rs\\d+\\s*$', SNPS)"

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
  # below is nonsense, so comment out
  # filter(`OR or BETA` > 0) %>%
  mutate(BETA = ifelse(!is.na(`OR or BETA`) & grepl('beta', `OR/RR/BETA`, ignore.case = TRUE), 
                        `OR or BETA`, log(`OR or BETA`))) %>%
  ##### khy
  # remove NaN which come from log of negative
  filter(!is.nan(BETA)) %>%
  
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
        ###########
        # what is meaning of gsub('^(rs\\d+)(\\D.*)?$', '\\1', snpname)?
        # [ref] https://stackoverflow.com/questions/37425019/gsub-only-part-of-pattern
        # [ans] "rs10001545:C:A" => "rs10001545"
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
  ###### khy
  # %dopar% -> %do%
  ###### khy
  # convert numeric
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
  mutate(SNP_OR = exp(BETA)) %>%
  (function(dat) {
    dat %>%
      (function(x) {
        ##### khy, debug
        fwrite(x, file = "before_merge_g1_samples.csv")
        #fwrite(dat, file = "before_checkLDPlink.rdata")
        cat("######## before merge g1_samples saved", sprintf('%s', Sys.time()), "\n")
        #### khy
        # "x %>% names" means to get column names of x object
        # below meaning: 
        # select data in column whose column name is matched with grepl pattern from x
        # "with = FALSE" give all data in column, "with = TRUE" give only column name
        #### khy
        # why below line is needed? there is no column name matched with '\\.(Initial|Replication)$' 
        # my guess: below 'reference/g1_samples.csv' file have column name matched with the pattern
        x[, !grepl('\\.(Initial|Replication)$', x %>% names), with = FALSE]
      }) %>%
      merge(fread('reference/g1_samples.csv', encoding = "UTF-8") %>% dplyr::select(-race, -Outcomes),
        by = c('PUBMEDID', 'LINK', 'code'), all.x = TRUE, all.y = FALSE)
  }) %>%
  mutate(
    East.Asian.Ancestry.Sample.Initial = 
      ifelse(East.Asian.Ancestry.Sample.Initial %in% c(0, NA),
        ifelse(!is.na(sample.EAS), sample.EAS, 0), East.Asian.Ancestry.Sample.Initial),
    East.Asian.Ancestry.Sample.Replication = 
      ifelse(East.Asian.Ancestry.Sample.Replication %in% c(0, NA),
        0, East.Asian.Ancestry.Sample.Replication)) %>%
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
  #### khy
  # study: (East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication > 0) * 15
  # operator precedence: (high) +, (low) >
  # pmax: pairwise max, pmin: pairwise min
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
    #### khy
    # "%dopar%" -> "%do%"
    #### khy
    #### MUST LEARN "foreach": 
    # "foreach(CODE = unique(dat$code), 
    #   .combine = function(...) rbindlist(list(...)), .multicombine=TRUE)"
    # [ref] https://stackoverflow.com/questions/17411223/r-foreach-with-combine-rbindlist
    foreach(CODE = unique(dat$code), 
            .combine = function(...) rbindlist(list(...)), .multicombine=TRUE) %do% {
      dat.code <- dat %>%
        filter(code == CODE)
      dat.code %>%
        arrange(desc(Point), desc(`P-VALUE`)) %>%
        mutate(beta_zero = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %do% {
          if (row == 1)
            return(FALSE)
          #### khy
          #### meaning of "dat.ld[SNP_A %in% snpname[1:(row - 1)], SNP_B]":
          # from dat.ld, select rows where "SNP_A" is in "snpname[1:(row - 1)]", and get "SNP_B" of those rows
          return(snpname[row] %in% dat.ld[SNP_A %in% snpname[1:(row - 1)], SNP_B])
        }) %>%
        mutate(BETA = ifelse(beta_zero, 0, BETA)) %>%
        mutate(SNP_OR = ifelse(beta_zero, 1, SNP_OR))
    }
  }) %>%
  group_by(code) %>%
  filter(!duplicated(snpname)) %>%
  mutate(order_id = 1:n()) %>%
  ungroup %>%
  mutate(gene_description = 'Not available') %>%
  mutate(title = 'NA') %>%
  rename(description = gene_description) %>%
  #### khy
  # "%dopar%" -> "%do%"
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
  rename(snp_or = SNP_OR) %>%
  mutate(snp_freq = pmin(snp_freq, 0.99999)) %>%
  mutate(snp_freq = pmax(snp_freq, 0.00001)) %>%
  mutate(aa_or = 1) %>%
  mutate(ab_or = ifelse(RECESSIVE, 1, snp_or)) %>%
  mutate(bb_or = ifelse(DOMINANT | RECESSIVE, snp_or, snp_or**2)) %>%
  mutate(lang = 'ko') %>%
  mutate(version = 1) %>%
  rename(ref_url = LINK) %>%
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

##### khy
# here, finally "unfiltered_meta" is obtained! long way
fwrite(unfiltered_meta, file = "unfiltered_meta.csv")
cat("##### unfiltered meta saved", sprintf("%s", Sys.time()), "\n")

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

#### khy, debug
fwrite(final_meta, file = "final_meta.csv")
cat("##### final meta(after bind_rows beta_zero) saved", sprintf("%s", Sys.time()), "\n")

final_meta %>%
  data.table %>%
  (function(dat) {
    setkey(dat, Bio.Serial)
    dat[data.table(Bio.Serial = 1:max(dat$Bio.Serial)) %>%
          (function(x) {
            setkey(x, Bio.Serial)
            x
          })
        ] %>%
        ######## khy, for debug
        ###### to know what is above mysterious code
        (function(dat) {
          fwrite(dat, file = "max_bio_serial.csv")
          cat("##### max bio serial saved", sprintf('%s', Sys.time()), "\n")
          dat
        }) %>%
      filter(!duplicated(Bio.Serial)) %>%
      mutate(OnMeta_DS = !is.na(code)) %>%
      mutate(`OnMeta_DS(del by LD)` = beta_zero) %>%
      dplyr::select(Bio.Serial, OnMeta_DS, `OnMeta_DS(del by LD)`, snpname) %>%
      arrange(Bio.Serial) %>%
      rename(OnMeta_snpname = snpname) %>%
      fwrite('G1_onmeta_joon.csv')
      cat("###### G1_onmeta_joon.csv saved", sprintf('%s', Sys.time()), "\n")
    dat
  }) %>%
  mutate(recessive = RECESSIVE %>% as.integer) %>%
  mutate(dominant = DOMINANT %>% as.integer) %>%
  mutate(`seq` = 1:n()) %>%
  dplyr::select(`seq`, code, order_id, title, snpname, gene, location, description,
                ref_allele, eff_allele, snp_freq, snp_or, snp_paf,
                aa, ab, bb, aa_or, ab_or, bb_or,
                lang, version, ref_url, TIMES, CHR, POS, ARM,
                position_percentage, recessive, dominant, beta_zero) %>%
  rename(chromosome = CHR, position = POS, arm = ARM) %>%
  mutate(gene = ifelse(is.na(gene), 'Intergenic', gene),
          location = ifelse(is.na(location), '-', location)) %>%
  mutate(ab_beta = log(ab_or),
          bb_beta = log(bb_or),
          ab_freq = 2 * snp_freq * (1 - snp_freq),
          bb_freq = snp_freq * snp_freq) %>%
  mutate(ab_beta_mean = ab_beta * ab_freq,
          bb_beta_mean = bb_beta * bb_freq,
          ab_beta_var = ab_beta * ab_beta * ab_freq * (1 - ab_freq),
          bb_beta_var = bb_beta * bb_beta * bb_freq * (1 - bb_freq)) %>%
  mutate(beta_mean = ab_beta_mean + bb_beta_mean,
          beta_var = ab_beta_var + bb_beta_var - 2 * ab_beta_mean * bb_beta_mean) %>%
  fwrite('G1_SNP_joon.csv')
  cat("###### G1_SNP_joon.csv saved", sprintf('%s', Sys.time()), "\n")

end_time <- Sys.time()
cat(sprintf("###### run time = %f min\n", 
    as.numeric(end_time - start_time, units = "mins")))
