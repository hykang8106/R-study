setwd("/Users/genoplan/metadata/meta_scripts/")

# /data/98.reference/11.meta_ref/Plink -> /Users/genoplan/metadata/reference/Plink
# /data/98.reference/09.dbsnp/ -> /Users/genoplan/metadata/reference/dbsnp/
# plink -> /Users/genoplan/tools/plink_mac_20210606/plink --bfile


knitr::opts_chunk$set(echo = TRUE)
options(encoding = "UTF-8")
Sys.setlocale("LC_ALL", "ko_KR.UTF-8")
Sys.setlocale("LC_CTYPE", "ko_KR.UTF-8")
# Sys.setenv(LANG = "ko_KR.UTF-8")
if (!require(data.table)) install.packages("data.table", repos = "https://Rdatatable.gitlab.io/data.table")
library(data.table)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(RSQLite)) install.packages('RSQLite')
library(RSQLite)
if (!require(doParallel)) install.packages('doParallel')
library(doParallel)
if (!require(openxlsx)) install.packages('openxlsx')
library(openxlsx)
if (!require(readr)) install.packages('readr')
library(readr)

load('cytogenic.band.length.RData')

#registerDoParallel(detectCores())
registerDoParallel(15)

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
  on.exit(unlink(sprintf('/Users/genoplan/metadata/reference/Plink/%s', inputFile)), add = TRUE)
  data.frame(snp = SNPList) %>%
    fwrite(sprintf('/Users/genoplan/metadata/reference/Plink/%s', inputFile), quote = FALSE, col.names = FALSE)
  outputFile <- sample(0:9, size = 20, replace = TRUE) %>% paste(collapse = '')
  on.exit(unlink(sprintf('/Users/genoplan/metadata/reference/Plink/%s.*', outputFile)), add = TRUE)
  system(sprintf('/Users/genoplan/tools/plink_mac_20210606/plink --bfile /Users/genoplan/metadata/reference/Plink/ref --memory 1536 --r2 --ld-window-kb 1000 --ld-window 5000 --ld-window-r2 0.7 --ld-snp-list /Users/genoplan/metadata/reference/Plink/%s --out /Users/genoplan/metadata/reference/Plink/%s', inputFile, outputFile), intern = TRUE)
  return(tryCatch(fread(sprintf('/Users/genoplan/metadata/reference/Plink/%s.ld', outputFile))[, .(SNP_A, SNP_B)], error = function(e) data.table(SNP_A = character(), SNP_B = character())))
}

# checkAlleleDbsnp
checkAlleleDbsnp <- function(SNP, Allele) {
  con <- dbConnect(SQLite(), '/Users/genoplan/metadata/reference/dbsnp/dbsnp.grch37.sqlite')
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

# findAltDbsnp
findAltDbsnp <- function(SNP, Allele) {
  if (grepl('^[ACGT]$', Allele, ignore.case = TRUE))
    return(Allele)
  con <- dbConnect(SQLite(), '/Users/genoplan/metadata/reference/dbsnp/dbsnp.grch37.sqlite')
  on.exit(dbDisconnect(con), add = TRUE)
  res <- dbGetQuery(con, sprintf('select ALT from dbsnp where `NAME`="%s" limit 1', SNP)) %>%
    unlist
  if (length(res) > 0)
    return(res)
  else
    return(Allele)
}

# 연속변수 평균값 처리
code.continuous <- fread('G1_Continuous_Mean.csv') %>%
  (function(dat) {
    dat[, .(code, PUBMEDID, LINK, mean, `DISEASE/TRAIT`)] %>%
      merge(dat %>%
              filter(!grepl('부적합', `Comment`)) %>%
              group_by(code) %>%
              summarize(n = sum(!is.na(mean))) %>%
              filter(n > 0), by = c('code'), all.x = FALSE, all.y = TRUE)
  }) %>%
  filter(!is.na(mean)) %>%
  mutate(LINK = ifelse(is.na(LINK), '', LINK)) %>%
  filter(!is.na(code)) %>%
  mutate(LINK = ifelse(!is.na(PUBMEDID) &
                         PUBMEDID != '' &
                         PUBMEDID != '0', sprintf('www.ncbi.nlm.nih.gov/pubmed/%s', PUBMEDID), LINK))

save(code.continuous, file = 'code.continuous.RData')



# Reference frequencey를 sqlite로 변환
unlink('/Users/genoplan/metadata/reference/Plink/ref.frq.sqlite')
fread('/Users/genoplan/metadata/reference/Plink/ref.frq') %>%
  (function(dat) {
    con <- dbConnect(SQLite(), '/Users/genoplan/metadata/reference/Plink/ref.frq.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    dat %>%
      dbWriteTable(conn = con, name = 'freq', overwrite = TRUE)
    dbExecute(conn = con, 'create index idx_freq_SNP_A1_A2 on freq (SNP, A1, A2)')
  })

# 바이오팀 데이터 가공
load('code.continuous.RData')
code.continuous.list <- fread('continuous.csv')$code %>%
  unique
unfiltered_meta <- fread('g1.csv') %>%
  # Bio.Serial column: Assign serial Number for each SNP
  mutate(Bio.Serial = 1:n()) %>% 
  filter(!is.na(code)) %>%
  # Checking history from BS : 1st check is '0'         
  filter(`1st check` == 'O') %>%
  
  # Substitute RSID col. with SNP col., if RSID col. is empty and SNP col. starts with 'rs'     
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
  
  # P-value to numeric (warning:NAs)
  mutate(`P-VALUE` = as.numeric(`P-VALUE`)) %>%
  
  # Select columns required to calculate
  dplyr::select(code, SNP, sample.EAS, sample.OTHER, sample.TOTAL,
                `OR/RR/BETA`, `OR or BETA`, Effect, `P-VALUE`, PUBMEDID, LINK, `Do/Re`, Bio.Serial) %>%
  
  # Rename column as readable in R
  rename('OR or BETA' = 'OR or BETA',
         'eff_allele' = 'Effect') %>%
  # Unify other form of effect of SNP to the word 'beta'; maintain 'OR'
  mutate(`OR/RR/BETA` = ifelse(!is.na(`OR/RR/BETA`) & grepl('(beta|% effect|coef|AUC|delta|mean|estimate|h2|u|z|ln\\(OR\\))', `OR/RR/BETA`, ignore.case = TRUE), 'BETA', `OR/RR/BETA`)) %>% #OR/RR/BETA의 정보를 beta 또는 OR으로 변경
  # OR or BETA to numeric (warnings:NAs)
  mutate(`OR or BETA` = as.numeric(`OR or BETA`)) %>%
  # OR to log(OR)
  mutate(BETA = ifelse(!is.na(`OR or BETA`) & grepl('beta', `OR/RR/BETA`, ignore.case = TRUE), `OR or BETA`, log(`OR or BETA`))) %>%
  
  # publication reference as pubmed link
  mutate(PUBMEDID = PUBMEDID %>% as.integer) %>%         #PUBMEDID를 integer로 변경
  mutate(LINK = ifelse(!is.na(PUBMEDID) &
                         PUBMEDID != '' &
                         PUBMEDID != '0', sprintf('www.ncbi.nlm.nih.gov/pubmed/%s', PUBMEDID), LINK)) %>% #LINK컬럼에 PUBMEDID가 있으면 pubmed link를 입력 아니면 그대로 LINK입력
  
  # ########################################################################
  mutate(DOMINANT = grepl('Do', `Do/Re`, ignore.case = TRUE)) %>%      #Do/Re컬럼에 "Do"이면 DOMINANT컬럼에 TRUE
  mutate(RECESSIVE = grepl('Re', `Do/Re`, ignore.case = TRUE)) %>%        #Do/Re컬럼에 "Re"이면 RECESSIVE컬럼에 TRUE
  mutate(eff_allele = toupper(eff_allele)) %>%                #eff_allele을 모두 대문자로 변경
  dplyr::select(-`Do/Re`) %>%                 #"Do/Re" 컬럼 제거
  mutate(eff_allele = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %dopar% findAltDbsnp(SNP[row], eff_allele[row])) %>%   #eff_allele컬럼에 dbsnp의 ref/alt 확인하여 없으면 strand 변경
  mutate(eff_allele = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %dopar% checkAlleleDbsnp(SNP[row], eff_allele[row])) %>%  #eff_allele컬럼에 dbsnp에서 없으면 alt sequence로 입력
  filter(eff_allele != '') %>%    #eff_allele이 없으면 제거
  data.table %>%    #data table 형태로 변경
  (function(dat) {
    setkey(dat, SNP)            #SNP를 key값으로 설정
    dat
  }) %>%
  (function(dat) {
    fread('/Users/genoplan/metadata/reference/Plink/ref.bim', header = FALSE)[, .(CHR = V1, snpname = V2, POS = V4, REF = V5, ALT = V6)] %>%  #reference bim file을 load
      (function(bim) {
        bim[, SNP := gsub('^(rs\\d+)(\\D.*)?$', '\\1', snpname)]                 #SNP 컬럼에 snp(rs number)만 뽑아서 입력
        setkey(bim, SNP)                                                       # SNP를 key값으로 설정
        bim[dat]
      })
  }) %>%
  filter(!is.na(CHR)) %>%                                         #Chromosome이 없으면 제거
  filter(REF==eff_allele | ALT==eff_allele) %>%                     #eff_allele이 REF 또는 ALT에 없으면 제거
  filter(!duplicated(Bio.Serial)) %>%                                 #Bio에서 부여한 number 중복 제거 (multi allele 때문에 생기는 현상)
  #filter(eff_allele == REF | eff_allele == ALT) %>%                  #!!!수정사항!!!! 순서변경
  mutate(ref_allele = ifelse(eff_allele == ALT, REF, ALT)) %>%          #ref_allele 컬럼에 eff_allele이 ALT와 같으면 REF를 입력 아니면 ALT를 입력
  mutate(aa = sprintf('%s%s', ref_allele, ref_allele)) %>%                #aa 컬럼에 "ref_allele""ref_allele"입력
  mutate(ab = ifelse(eff_allele > ref_allele,                               #ab 컬럼에 "ref_allele""alt_allele" 또는 "alt_allele""ref_allele" 입력
                     sprintf('%s%s', ref_allele, eff_allele),
                     sprintf('%s%s', eff_allele, ref_allele))) %>%
  mutate(bb = sprintf('%s%s', eff_allele, eff_allele)) %>%                     #bb 컬럼에 "alt_allele""alt_allele"입력
  ungroup %>%
  data.table %>%
  (function(dat) {
    registerDoParallel(detectCores())                                                    #현재 사용 가능한 Core 적용
    dat[, location := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %dopar% {            #mutiThread 사용하여 라인마다 처리
      if (!CHR[row] %in% unique(cytogenic.band.length$chr))                                   #cytogenic bandDB에 CHR이 존재하지 않으면 "-"
        return('-')                                                                             
      row.band <- cytogenic.band.length %>%                                                 #cytogenic bandDB에 CHR이 존재하면
        filter(chr == CHR[row]) %>%                                                         #cytogenic bandDB의 chr과 CHR이 같은것 선택
        filter(bp.start <= POS[row]) %>%                                                   #cytogenic bandDB의 bp.start 보다 큰 것 선택
        filter(bp.stop >= POS[row])                                                        #cytogenic bandDB의 bp.end 보다 작은 것 선택 하여 raw.band에 저장
      if (nrow(row.band) == 0)                                                            #row.band에 들어갈 데이터가 없으면
        return('-')                                                                       # "-"를 입력
      return(sprintf('%s%s%s',                                                              #chr(p/q)band return
                     CHR[row],
                     row.band$arm[1],
                     row.band$band[1]))
    }]
    dat[, ARM := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %dopar% {                
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
    dat[, position_percentage := foreach(row = 1:.N, .combine = c, .multicombine = TRUE) %dopar% {
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
    dat
  })  %>%
  mutate(gene = foreach(snp = SNP, .combine = c, .multicombine = TRUE) %dopar% {
    con <- dbConnect(SQLite(), '/Users/genoplan/metadata/reference/dbsnp/dbsnp.gene.grch37.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    res <- dbGetQuery(con, sprintf('select GENE from gene where SNP="%s"', snp))
    if (nrow(res) == 0)
      return('intergenic')
    else
      return(res$GENE[1])
  }) %>%
  mutate(LINK = ifelse(is.na(LINK), '', LINK)) %>%
  mutate(TIMES = foreach(row = 1:n(), .combine = c) %dopar% {
    return(!is.na(BETA[row]) &
             (!code[row] %in% code.continuous.list |
                (sum(code.continuous$code == code[row] &
                       code.continuous$LINK == LINK[row] &
                       !(is.na(code.continuous$mean) | code.continuous$mean == 0)) > 0)
             ))
  }) %>%
  mutate(BETA = foreach(row = 1:n(), .combine = c) %dopar% {
    if (sum(TIMES[code == code[row]]) == 0)
      return(0.001 * BETA[row] / abs(BETA[row]))
    else if (!code[row] %in% code.continuous.list)
      return(BETA[row])
    mean <- code.continuous$mean[code.continuous$code == code[row] & code.continuous$LINK == LINK[row]]
    mean <- mean %>% setdiff(NA)
    if (length(mean) == 0)
      return(0.001 * BETA[row] / abs(BETA[row]))
    
    mean <- mean[1]
    if (mean == 0)
      return(0.001 * BETA[row] / abs(BETA[row]))
    else if (mean < 0)
      return(BETA[row])
    else
      return(log((BETA[row] + mean) / mean))
  }) %>%
  filter(!is.na(BETA)) %>%
(function(dat) {
  dat %>%
    group_by(code) %>%
    summarize(TIMES = sum(TIMES, na.rm = TRUE)) %>%
    filter(TIMES == 0) %>%
    dplyr::select(code) %>%
    fwrite('G1_TIMES.csv')
  dat
}) %>%
  mutate(SNP_OR = exp(BETA)) %>%
  (function(dat) {
    dat %>%
      (function(x) {
        x[,!grepl('\\.(Initial|Replication)$', x %>% names), with = FALSE] 
      }) %>%
      merge(fread('g1_samples.csv') %>% dplyr::select(-race, -Outcomes), by = c('PUBMEDID', 'LINK', 'code'), all.x = TRUE, all.y = FALSE)
  }) %>%
  mutate(East.Asian.Ancestry.Sample.Initial = ifelse(East.Asian.Ancestry.Sample.Initial %in% c(0, NA),
                                                     ifelse(!is.na(sample.EAS), sample.EAS, 0),
                                                     East.Asian.Ancestry.Sample.Initial),
         East.Asian.Ancestry.Sample.Replication = ifelse(East.Asian.Ancestry.Sample.Replication %in% c(0, NA),
                                                         0,
                                                         East.Asian.Ancestry.Sample.Replication)) %>%
  mutate(Asian.Ancestry.Sample.Initial = ifelse(Asian.Ancestry.Sample.Initial %in% c(0, NA),
                                                ifelse(is.na(East.Asian.Ancestry.Sample.Initial), 0, East.Asian.Ancestry.Sample.Initial),
                                                Asian.Ancestry.Sample.Initial),
         Asian.Ancestry.Sample.Replication = ifelse(Asian.Ancestry.Sample.Replication %in% c(0, NA),
                                                    ifelse(is.na(East.Asian.Ancestry.Sample.Replication), 0, East.Asian.Ancestry.Sample.Replication),
                                                    Asian.Ancestry.Sample.Replication),
         European.Ancestry.Sample.Initial = ifelse(European.Ancestry.Sample.Initial %in% c(0, NA),
                                                   0,
                                                   European.Ancestry.Sample.Initial),
         European.Ancestry.Sample.Replication = ifelse(European.Ancestry.Sample.Replication %in% c(0, NA),
                                                       0,
                                                       European.Ancestry.Sample.Replication)) %>%
  mutate(Total.Ancestry.Sample.Initial = ifelse(Total.Ancestry.Sample.Initial %in% c(0, NA),
                                                ifelse(!is.na(sample.TOTAL), sample.TOTAL, Asian.Ancestry.Sample.Initial + European.Ancestry.Sample.Initial),
                                                Total.Ancestry.Sample.Initial),
         Total.Ancestry.Sample.Replication = ifelse(Total.Ancestry.Sample.Replication %in% c(0, NA),
                                                    Asian.Ancestry.Sample.Replication + European.Ancestry.Sample.Replication,
                                                    Total.Ancestry.Sample.Replication)) %>%
  mutate(Point = (ifelse(is.na(East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication), 0, (East.Asian.Ancestry.Sample.Initial + East.Asian.Ancestry.Sample.Replication > 0) * 15) +
                    ifelse(is.na(Asian.Ancestry.Sample.Initial + Asian.Ancestry.Sample.Replication), 0, (Asian.Ancestry.Sample.Initial + Asian.Ancestry.Sample.Replication > 0) * 30) +
                    ifelse(is.na(European.Ancestry.Sample.Initial + European.Ancestry.Sample.Replication), 0, (European.Ancestry.Sample.Initial + European.Ancestry.Sample.Replication > 0) * 10) +
                    ifelse(is.na(Total.Ancestry.Sample.Initial + Total.Ancestry.Sample.Replication), 0, pmax(pmin(log10(Total.Ancestry.Sample.Initial + Total.Ancestry.Sample.Replication), 5), 0) * 10) +
                    ifelse(is.na(`P-VALUE`), 0, pmin(-log10(`P-VALUE`), 10))) / 1.15
  ) %>%
  (function(dat) {
    dat.ld <- checkLDPlink(dat$snpname %>% unique %>% setdiff(c(NA, '')))
    foreach(CODE=unique(dat$code), .combine=function(...) rbindlist(list(...)), .multicombine=TRUE) %dopar% {
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
  group_by(code) %>%
  filter(!duplicated(snpname)) %>%
  mutate(order_id = 1:n()) %>%
  ungroup %>%
  mutate(gene_description = 'Not available') %>%
  mutate(title = 'NA') %>%
  rename(description = gene_description) %>%
  mutate(snp_freq = foreach(row = 1:n(), .combine = c, .multicombine = TRUE) %dopar% {
    con <- dbConnect(SQLite(), '/Users/genoplan/metadata/reference/Plink/ref.frq.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    res <- dbGetQuery(con, sprintf('select * from freq where SNP="%s" and ((A1="%s" and A2="%s") or (A1="%s" and A2="%s"))',
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
  mutate(snp_paf = jointPAF(cbind(ifelse(ab_or > 1, 2 * snp_freq * (1 - snp_freq) * (ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (ab_or - 1)), 2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1) / (1 + 2 * snp_freq * (1 - snp_freq) * (1 / ab_or - 1))),
                                  ifelse(bb_or > 1, snp_freq * snp_freq * (bb_or - 1) / (1 + snp_freq * snp_freq * (bb_or - 1)), (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1) / (1 + (1 - snp_freq) * (1 - snp_freq) * (1 / bb_or - 1)))))) #%>%


beta_zero <- subset(unfiltered_meta, unfiltered_meta$beta_zero==TRUE)

final_meta  <- subset(unfiltered_meta, unfiltered_meta$beta_zero==FALSE) %>%
  (function(dat) {
    foreach(CODE = unique(dat$code), .combine = function(...) rbindlist(list(...)), .multicombine = F) %dopar% {
      for (cutoff in seq(0.3, 1, 0.02)) {
        res <- dat %>%
          filter(code == CODE) %>%
          filter(snp_paf <= cutoff )
        if (nrow(res) > 1)
          break
      }
      return(res)
    }
  }) %>% bind_rows(beta_zero)

final_meta %>%
  data.table %>%
  (function(dat) {
    setkey(dat, Bio.Serial)
    dat[data.table(Bio.Serial = 1:max(dat$Bio.Serial)) %>%
          (function(x) {
            setkey(x, Bio.Serial)
            x
          })] %>%
      filter(!duplicated(Bio.Serial)) %>%
      mutate(OnMeta_DS = !is.na(code)) %>%
      mutate(`OnMeta_DS(del by LD)` = beta_zero) %>%
      dplyr::select(Bio.Serial, OnMeta_DS, `OnMeta_DS(del by LD)`, snpname) %>%
      arrange(Bio.Serial) %>%
      rename(OnMeta_snpname = snpname) %>%
      fwrite('G1_onmeta_joon.csv')
    dat
  }) %>%
  mutate(recessive = RECESSIVE %>% as.integer) %>%
  mutate(dominant = DOMINANT %>% as.integer) %>%
  mutate(`seq` = 1:n()) %>%
  dplyr::select(`seq`, code, order_id, title, snpname, gene, location, description, ref_allele, eff_allele, snp_freq, snp_or, snp_paf, aa, ab, bb, aa_or, ab_or, bb_or, lang, version, ref_url, TIMES, CHR, POS, ARM, position_percentage, recessive, dominant, beta_zero) %>%
  rename(chromosome = CHR, position = POS, arm = ARM) %>%
  mutate(gene = ifelse(is.na(gene), 'Intergenic', gene), location = ifelse(is.na(location), '-', location)) %>%
  mutate(ab_beta = log(ab_or), bb_beta = log(bb_or), ab_freq = 2 * snp_freq * (1 - snp_freq), bb_freq = snp_freq * snp_freq) %>%
  mutate(ab_beta_mean = ab_beta * ab_freq, bb_beta_mean = bb_beta * bb_freq, ab_beta_var = ab_beta * ab_beta * ab_freq * (1 - ab_freq), bb_beta_var = bb_beta * bb_beta * bb_freq * (1 - bb_freq)) %>%
  mutate(beta_mean = ab_beta_mean + bb_beta_mean, beta_var = ab_beta_var + bb_beta_var - 2 * ab_beta_mean * bb_beta_mean) %>%
  fwrite('G1_SNP_joon.csv')


# test_meta <- final_meta
# 
# test_meta <- test_meta %>% mutate(filter_key = paste(code, gene, location, sep = "@"))
# filtered_key <- unique(test_meta$filter_key)
# 
# #for (i in 1:length(filter_key)){
# 
# tmp <- NULL
# tmp_table <- NULL
# for (i in 1:length(filtered_key)){
#   tmp[[i]] <- subset(test_meta, test_meta$filter_key == filtered_key[i]) %>% arrange(beta_zero, desc(position))
#   if(nrow(subset(tmp[[i]], tmp[[i]]$beta_zero == FALSE)) > 2){
#     tmp_table[[i]] <- subset(tmp[[i]], tmp[[i]]$beta_zero == FALSE)
#   }else{
#     tmp_table[[i]] <- tmp[[i]][c(1:3),] %>% dplyr::filter(!is.na(code))
#   } 
# }
# 
# zz <- ldply(tmp_table, rbind)
# 
# 
