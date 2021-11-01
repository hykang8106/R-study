
'''
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

# "chromosome" is same as "cytogenic_band_length['chr'].unique()"
# to reduce run time, may use "chromosome"
chromosome = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', \
  '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

def location_from_cytogenic_band(cytogenic_band_length, CHR, POS):

  # global chromosome

  # type(cytogenic_band_length['chr'][0] is 'str'
  # type(CHR) is 'int64'
  # type(POS) is 'int64'

  # if str(CHR) not in cytogenic_band_length['chr'].unique():
  if str(CHR) not in chromosome:
    return '-'

  cyto_filter = (cytogenic_band_length['chr'] == str(CHR))
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  cyto_filter = (cytogenic_band_length['bp.start'] <= POS)
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  cyto_filter = (cytogenic_band_length['bp.stop'] >= POS)
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  if len(cytogenic_band_length) == 0:
    return '-'
  
  return (str(CHR) + cytogenic_band_length['arm'][0] + str(cytogenic_band_length['band'][0]))


def ARM_from_cytogenic_band(cytogenic_band_length, CHR, POS):

  # global chromosome

  if str(CHR) not in chromosome:
    return '-'

  cyto_filter = (cytogenic_band_length['chr'] == str(CHR))
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  cyto_filter = (cytogenic_band_length['bp.start'] <= POS)
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  cyto_filter = (cytogenic_band_length['bp.stop'] >= POS)
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  if len(cytogenic_band_length) == 0:
    return '-'
  
  return cytogenic_band_length['arm'][0]


def position_percentage_from_cytogenic_band(cytogenic_band_length, CHR, ARM, POS):

  # global chromosome

  if str(CHR) not in chromosome:
    return 0

  cyto_filter = (cytogenic_band_length['chr'] == str(CHR))
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  cyto_filter = (cytogenic_band_length['arm'] == ARM)
  cytogenic_band_length = cytogenic_band_length[cyto_filter].reset_index(drop=True)

  if len(cytogenic_band_length) == 0:
    return 0

  return (POS - min(cytogenic_band_length['bp.start'])) / \
    (max(cytogenic_band_length['bp.stop']) - min(cytogenic_band_length['bp.start']))

