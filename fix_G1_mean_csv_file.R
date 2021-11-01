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
#writeLines(new_text_line, 'fixed_G1_Continuous_Mean.csv')
#write.csv(new_text_line, 'test.csv', row.names = FALSE, na = "")
writeLines(new_text_line, 'reference/fixed_G1_Continuous_Mean.csv')
cat("##### fixed G1_Continuous_Mean csv file was created\n")
