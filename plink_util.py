# plink utility

import pandas as pd
import numpy as np
import sqlite3 as sql
import os

from file_db_define import *

def write_plink_ref_freq(plink_ref_freq_file, plink_ref_freq_db):

   # separation char is one more space, so sep='\s+', and grep expression need 'python' engine
   plink_ref_freq = \
        pd.read_csv(plink_ref_freq_file, encoding="UTF-8", sep='\s+', engine='python')
   #ref_freq = pd.read_csv(plink_ref_freq_file, encoding="UTF-8", sep='\t')
   print(plink_ref_freq.shape)
   print(plink_ref_freq.head())

   # to have database in ram, use below ":memory:"
   # con = sql.connect(":memory:")

   # remove database if exist
   if os.path.isfile(plink_ref_freq_db):
      os.remove(plink_ref_freq_db)

   con = sql.connect(plink_ref_freq_db)
   cur = con.cursor()

   # below line is not needed because database was removed already
   # if database is not removed in advance, use below line
   # cur.execute("DROP TABLE IF EXISTS freq")

   print("#### write plink_ref_freq into freq table\n")

   # data type: [ref] https://docs.python.org/ko/3/library/sqlite3.html
   cur.execute("CREATE TABLE freq(CHR INTEGER, SNP TEXT, A1 TEXT, A2 TEXT, MAF REAL, NCHROBS INTEGER)")
   # cur.execute("CREATE TABLE freq(CHR tinyint, SNP text, A1 text, A2 text, MAF real, NCHROBS smallint)")

   # "executemany" need "list of tuple" input, so convert dataframe to list of tuple": 
   # [ref] https://www.kite.com/python/answers/how-to-convert-a-pandas-dataframe-into-a-list-of-tuples-in-python
   cur.executemany("INSERT INTO freq VALUES (?, ?, ?, ?, ?, ?)", \
      list(plink_ref_freq.to_records(index=False)))

   cur.execute("CREATE INDEX idx_freq_SNP_A1_A2 ON freq(SNP, A1, A2)")

   con.commit()
   con.close()

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

def query_plink_ref_freq(plink_ref_freq_db, snpname, ref_allele, eff_allele):

   con = sql.connect(plink_ref_freq_db)
   cur = con.cursor()
   query_string = \
        'SELECT * FROM freq WHERE SNP="{}" and ((A1="{}" and A2="{}") or (A1="{}" and A2="{}"))'.\
        format(snpname, ref_allele, eff_allele, eff_allele, ref_allele)
   cur.execute(query_string)
   rows = cur.fetchall()
   # "rows" type = list of tuple
   # "rows" example:
   # CHR                                   SNP               A1   A2  MAF         NCHROBS
   # [(b'\x01\x00\x00\x00\x00\x00\x00\x00', 'rs9645429:G:A', 'A', 'G', 0.001823, b'\xb8\x19\x00\x00\x00\x00\x00\x00')]
   con.close()

   if len(rows) == 0:
      return np.NaN
   else:
      rows = rows[0]
      return np.where(rows[2] == eff_allele, rows[4], 1 - rows[4])

def test_query_plink_ref_freq(plink_ref_freq_db):

   print("### query plink_ref_freq_db\n")

   con = sql.connect(plink_ref_freq_db)
   cur = con.cursor()
   query_string = \
        'SELECT * FROM freq WHERE SNP="{}" and ((A1="{}" and A2="{}") or (A1="{}" and A2="{}"))'.\
        format("rs9645429:G:A", "A", "G", "A", "G")
   cur.execute(query_string)
   rows = cur.fetchall() 
   con.close()
   for row in rows: 
        print(row)
        # print example:
        # "CHR", "NCHROBS" column(integer type) show byte
        # CHR                                   SNP               A1   A2  MAF         NCHROBS
        # (b'\x01\x00\x00\x00\x00\x00\x00\x00', 'rs9645429:G:A', 'A', 'G', 0.001823, b'\xb8\x19\x00\x00\x00\x00\x00\x00')

   return rows


if __name__ == "__main__":

   test_write_plink_ref_freq = False

   if test_write_plink_ref_freq:
      write_plink_ref_freq(plink_ref_freq_file, plink_ref_freq_db)

   rows = test_query_plink_ref_freq(plink_ref_freq_db)