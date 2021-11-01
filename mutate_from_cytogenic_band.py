


############ check: what "ungroup" in R code do?

def mutate_from_cytogenic_band(bim, cytogenic_band_length_file):

   import pandas as pd

   cytogenic_band_length = pd.read_csv(cytogenic_band_length_file, encoding="UTF-8")
   print(cytogenic_band_length.shape)
   print(cytogenic_band_length.head())

   '''
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
   '''

   # sys.exit()

   for i in range(len(bim)):

      bim.loc[i, 'location'] = \
         location_from_cytogenic_band(cytogenic_band_length, bim['CHR'][i], bim['POS'][i])

      bim.loc[i, 'ARM'] = \
         ARM_from_cytogenic_band(cytogenic_band_length, bim['CHR'][i], bim['POS'][i])

      bim.loc[i, 'position_percentage'] = \
         position_percentage_from_cytogenic_band(cytogenic_band_length, \
         bim['CHR'][i], bim['ARM'][i], bim['POS'][i])

   bim.to_csv('cytogenic_refered_py.csv', encoding='utf-8', index=False)
