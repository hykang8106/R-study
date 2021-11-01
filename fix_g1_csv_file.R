#### khy
# fix "g1.csv":
# (1) if any, remove double quote in "Korean" column
# (2) save into new file

setwd("C:/Users/hykang/Desktop/R_study")
cat("## current working directory =", getwd(), "\n")

# for "fread", "fwrite"
if (!require(data.table)) install.packages("data.table")
library(data.table)

# for "mutate", "%>%(pipeline)"
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

dat <- fread('reference/g1.csv', encoding = 'UTF-8') %>% 
  #mutate(code = gsub('A', 'Z', code))
  # remove double quote character in "Korean" column
  # this make mis-aligned data in "SNPS", "RSID" to be shifted into right column
  # for mis-aligned example, see line 2098 in original "g1.csv" file
  mutate(Korean = gsub("\"", "", Korean))

write.csv(dat, 'reference/fixed_g1.csv', row.names = FALSE, na = "")

cat("##### fixed g1 csv file was created\n")
