# to write python version of "meta_test_joon_win10.R"
# test sql
# [ref] https://docs.python.org/3/library/sqlite3.html

import pandas as pd
import sqlite3 as sql
import sys, os

plink_ref_freq_file = "C:/Users/hykang/Desktop/R_study/reference/Plink/ref.frq"
# "test.ref.frq.sqlite": not overwrite original db, "ref.frq.sqlite"
plink_ref_freq_db = "C:/Users/hykang/Desktop/R_study/reference/Plink/test.ref.frq.sqlite"

dbsnp_grch37_db = "C:/Users/hykang/Desktop/R_study/reference/dbsnp/dbsnp.grch37.sqlite"
dbsnp_gene_grch37_db = "C:/Users/hykang/Desktop/R_study/reference/dbsnp/dbsnp.gene.grch37.sqlite"

# plink test take long time, so to skip, set to "False"
test_plink = False
# test_plink = True

if test_plink:

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

    print("### query plink_ref_freq_db\n")

    con = sql.connect(plink_ref_freq_db)
    cur = con.cursor()
    query_string = \
        'SELECT * FROM freq WHERE SNP="{}" and ((A1="{}" and A2="{}") or (A1="{}" and A2="{}"))'.\
        format("rs9645429:G:A", "A", "G", "A", "G")
    cur.execute(query_string)
    rows = cur.fetchall() 
    for row in rows: 
        print(row)
        # print example:
        # "CHR", "NCHROBS" column(integer type) show byte
        # (b'\x01\x00\x00\x00\x00\x00\x00\x00', 'rs9645429:G:A', 'A', 'G', 0.001823, b'\xb8\x19\x00\x00\x00\x00\x00\x00')

    con.close()

# to test query in "findAltDbsnp" function
# res <- dbGetQuery(con, sprintf('select ALT from dbsnp where `NAME`="%s" limit 1', SNP))
print("### query dbsnp_grch37_db to get ALT\n")

snp = 'rs10099338'
con = sql.connect(dbsnp_grch37_db)
cur = con.cursor()
query_string = 'SELECT ALT FROM dbsnp WHERE NAME="{}" LIMIT 1'.format(snp)
cur.execute(query_string)
rows = cur.fetchall() 
for row in rows: 
    print(row)
    # print example:
    # ('G',)

con.close()

# to test query in "checkAlleleDbsnp" function
# res <- dbGetQuery(con, sprintf('select REF, ALT from dbsnp where `NAME`="%s"', SNP))
print("### query dbsnp_grch37_db to get REF, ALT\n")

snp = 'rs10032216'
con = sql.connect(dbsnp_grch37_db)
cur = con.cursor()
query_string = 'SELECT REF, ALT FROM dbsnp WHERE NAME="{}"'.format(snp)
cur.execute(query_string)
rows = cur.fetchall() 
for row in rows: 
    print(row)
    # print example:
    # ('T', 'A')
    # ('T', 'C')

con.close()

# to test query for "GENE"
# res <- dbGetQuery(con, sprintf('select GENE from gene where SNP="%s"', snp))
print("### query dbsnp_gene_grch37_db to get GENE\n")

con = sql.connect(dbsnp_gene_grch37_db)
cur = con.cursor()

# below snp give "('SREK1IP1',)" result
snp = "rs10075967"
# below snp give "empty" result
# snp = "rs10032216"
query_string = 'SELECT GENE FROM gene WHERE SNP="{}"'.format(snp)
cur.execute(query_string)
rows = cur.fetchall() 
if len(rows):
    for row in rows: 
        print(row)
else:
    print("#### GENE query result is empty\n")

con.close()

# sys.exit()

'''
con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.gene.grch37.sqlite')
on.exit(dbDisconnect(con), add = TRUE)
res <- dbGetQuery(con, sprintf('select GENE from gene where SNP="%s"', snp))

con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.grch37.sqlite')
on.exit(dbDisconnect(con), add = TRUE)
res <- dbGetQuery(con, sprintf('select ALT from dbsnp where `NAME`="%s" limit 1', SNP))

con <- dbConnect(SQLite(), 'reference/dbsnp/dbsnp.grch37.sqlite')
on.exit(dbDisconnect(con), add = TRUE)
res <- dbGetQuery(con, sprintf('select REF, ALT from dbsnp where `NAME`="%s"', SNP))
'''

