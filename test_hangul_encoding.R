# test hangul encoding

if (!require(readr)) install.packages('readr')
library(readr)

if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

if (!require(data.table)) install.packages("data.table")
library(data.table)

Sys.getlocale()

file_name <- "C:/Users/hykang/Desktop/R_study/reference/fixed_G1_Continuous_Mean.csv"

# [ref] http://drtagkim.blogspot.com/2019/02/r-csv.html
my_data <- read_csv(file_name, locale = locale("ko", encoding = "euc-kr"),
                    show_col_types = "FALSE")
my_data <- data.table(my_data)
names(my_data)
head(my_data)
dim(my_data)

f <- filter(my_data, my_data$code == "AAA")
##### why hangul text not working?
#f <- dplyr::filter(my_data, my_data$`Comment` == "결과 부적합")
head(f)

hh <- my_data$Comment
class(hh)
str(hh)
#dim(hh)

quit(save = "no")

#### khy
# fix "G1_Continuous_Mean.csv":
# (1) if any, remove double quote in line
# (2) save into new file

setwd("C:/Users/hykang/Desktop/R_study")
cat("## current working directory =", getwd(), "\n")

# read line from G1 mean csv file
#text_line <- readLines('reference/G1_Continuous_Mean.csv', encoding='EUC-KR')
text_line <- readLines('reference/G1_Continuous_Mean.csv', encoding='UTF-8')
#text_line <- readLines('reference/G1_Continuous_Mean.csv')
line_length <- length(text_line)

# allocate character vector
new_text_line <- vector('character', line_length)
for(i in 1:line_length) {
  # remove double quote
  new_text_line[i] <- gsub('\"', '', text_line[i])
}

# write into fixed G1 mean csv file
#write.table(new_text_line, 'test.csv', row.names = FALSE, na = "", sep = ",")
writeLines(new_text_line, 'test.csv')
cat("##### fixed G1_Continuous_Mean csv file was created\n")

quit(save = "no")

#################################################
######## rewritten(210909), NOT working, why?
#### try remove "," in "Korean" column
# "fix_g1_csv_file.R" work

# for "fread", "fwrite"
if (!require(data.table)) install.packages("data.table")
library(data.table)

# for "mutate", "%>%(pipeline)"
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

dat <- fread('reference/G1_Continuous_Mean.csv', encoding = 'UTF-8') %>% 
  #mutate(code = gsub('A', 'Z', code))
  # remove double quote character in "Korean" column
  # this make mis-aligned data in "SNPS", "RSID" to be shifted into right column
  # for mis-aligned example, see line 2098 in original "g1.csv" file
  mutate(Korean = gsub("\"", "", Korean))

write.csv(dat, 'fixed_G1_Continuous_Mean.csv', row.names = FALSE, na = "")

cat("##### fixed G1_Continuous_Mean csv file was created\n")
