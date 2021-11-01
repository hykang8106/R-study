# dbsnp query function: "findAltDbsnp", "checkAlleleDbsnp", "findGeneDbsnp"
# "findGeneDbsnp" function is new in python code (not in R code)

import pandas as pd
import sqlite3 as sql
import re, sys, os

dbsnp_grch37_db = "C:/Users/hykang/Desktop/R_study/reference/dbsnp/dbsnp.grch37.sqlite"
dbsnp_gene_grch37_db = "C:/Users/hykang/Desktop/R_study/reference/dbsnp/dbsnp.gene.grch37.sqlite"

'''
# rev.strand
rev.strand <- function(x) {
  foreach(x1 = strsplit(x, '')[[1]], .combine = c) %do% {
    list(A = 'T', T = 'A', C = 'G', G = 'C', D = 'D', I = 'I', '-' = '-')[[x1]]
  } %>%
    paste(collapse = '')
}
'''

def rev_strand(x):

  # R code used "named list": list(A = 'T', T = 'A', C = 'G', G = 'C', D = 'D', I = 'I', '-' = '-')
  # python code use "dict" (key = x, value = rev_x)
  dict_rev = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'D':'D', 'I':'I', '-':'-'}

  result = ''
  for x1 in list(x):
    # print(x1)
    x1_rev = dict_rev[x1]
    # print(x1_rev)
    result += x1_rev

  return result

'''
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
'''

def findAltDbsnp(SNP, Allele):

  # almost case return with "Allele": 31808 of 31814
  # database query is very rare: only 6 of 31814
  # this is right?

  if (re.search('^[ACGT]$', str(Allele), re.IGNORECASE) is not None):
    return Allele
    
  con = sql.connect(dbsnp_grch37_db)
  cur = con.cursor()
  query_string = 'SELECT ALT FROM dbsnp WHERE NAME="{}" LIMIT 1'.format(SNP)
  cur.execute(query_string)
  rows = cur.fetchall()
  con.close()

  # print(rows, type(rows))
  # print example: [('G',)], type = list of tuple
  if len(rows):
    # 1st element of list, 1st element of tuple: too naive?
    res = rows[0][0]
    return res
  else:
    return Allele

'''
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
'''

def checkAlleleDbsnp(SNP, Allele):

  con = sql.connect(dbsnp_grch37_db)
  cur = con.cursor()

  query_string = 'SELECT REF, ALT FROM dbsnp WHERE NAME="{}"'.format(SNP)
  cur.execute(query_string)
  rows = cur.fetchall() 
  con.close()

  # print(rows, type(rows))
  # print example: [('T', 'A'), ('T', 'C')], type = list of tuple

  # nested list comprehension: equivalent to "unlist" in R code
  res = [item for t in rows for item in t]
  # print(res, type(res))
  # print example = ['T', 'A', 'T', 'C']

  if Allele in res:
    return Allele
  elif rev_strand(Allele) in res:
    # print('### "rev_strand" function is called')
    return rev_strand(Allele)
  else:
    return ''

'''
con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.gene.grch37.sqlite')
    on.exit(dbDisconnect(con), add = TRUE)
    res <- dbGetQuery(con, sprintf('select GENE from gene where SNP="%s"', snp))
    if (nrow(res) == 0)
      return('intergenic')
    else
      return(res$GENE[1])
'''

def findGeneDbsnp(snp):

  con = sql.connect(dbsnp_gene_grch37_db)
  cur = con.cursor()

  # below snp give "('SREK1IP1',)" result
  # snp = "rs10075967"
  # below snp give "empty" result
  # snp = "rs10032216"
  query_string = 'SELECT GENE FROM gene WHERE SNP="{}"'.format(snp)
  cur.execute(query_string)
  rows = cur.fetchall() 
  con.close()

  # print(rows)
  if len(rows):
    res = rows[0][0]
    return res
  else:
    return 'intergenic'

# below "if" enable to test function before import:
# to test function, "python query_dbsnp.py" in windows command prompt
if __name__ == "__main__":

  snp = 'rs10099338'
  # for test, give invalid 'K'
  allele = 'K'
  # allele: one of 'A','C','G','T'
  # allele = 'A'
  res = findAltDbsnp(snp, allele)
  print('### "findAltDbsnp" return with', res, '\n')

  snp = 'rs10032216'
  # allele: one of 'A','C','G','T'
  # allele = 'T'
  # to test "rev_strand" function
  allele = 'G'
  res = checkAlleleDbsnp(snp, allele)
  print('#### "checkAlleleDbsnp" return with', res, '\n')

  x = 'ATCGDI-'
  # x = 'A'
  rev_x = rev_strand(x)
  print('#### "rev_strand":', x, '=>', rev_x, '\n')

  # below snp give "('SREK1IP1',)" result
  snp = "rs10075967"
  # below snp give "empty" result
  # snp = "rs10032216"
  res = findGeneDbsnp(snp)
  print('#### "findGeneDbsnp" return with', res, '\n')

