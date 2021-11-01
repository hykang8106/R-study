##### main routine
##### program run time(hh:mm:ss) = 00:37:13

import time

from get_code_continuous_mean import *
from plink_util import *
from process_g1 import *
from merge_ref_bim import *
from process_bim import *
from merge_g1_samples import *
from get_final_meta import *

from file_db_define import *

start_time = time.time()

code_continuous = get_code_continuous_mean(g1_continuous_mean_file)

write_plink_ref_freq(plink_ref_freq_file, plink_ref_freq_db)

g1 = process_g1(g1_file)

bim = merge_ref_bim(ref_bim_file, g1)

bim = mutate_from_cytogenic_band(bim, cytogenic_band_length_file)

bim = process_bim(bim, continuous_file, code_continuous)

unfiltered_meta = merge_g1_samples(bim, g1_samples_file)

final_meta, beta_zero = get_final_meta(unfiltered_meta)

end_time = time.time()

print('### program run time(hh:mm:ss) =', \
time.strftime('%H:%M:%S', time.gmtime(end_time - start_time)))

########################################
#### below is not written

# get_G1_onmeta()

# get_G1_SNP()

