# file and database define

# change below as user need
work_dir = "C:/Users/hykang/Desktop/R_study/"

g1_continuous_mean_file = work_dir + "reference/fixed_G1_Continuous_Mean.csv"
code_continous_file = work_dir + "code_continuous_py.csv"

plink_ref_freq_file = work_dir + "reference/Plink/ref.frq"
# "test.ref.frq.sqlite": not overwrite original db, "ref.frq.sqlite"
plink_ref_freq_db = work_dir + "reference/Plink/test.ref.frq.sqlite"

g1_file = work_dir + "reference/fixed_g1.csv"
before_bim_file = work_dir + "beforebim_py.csv"

continuous_file = work_dir + "reference/continuous.csv"
ref_bim_file = work_dir + "reference/Plink/ref.bim"
cytogenic_band_length_file = work_dir + "reference/cytogenic.band.length.csv"

G1_TIMES_file = work_dir + "G1_TIMES_py.csv"
before_merge_g1_samples_file = work_dir + "before_merge_g1_samples_py.csv"

dbsnp_grch37_db = work_dir + "reference/dbsnp/dbsnp.grch37.sqlite"
dbsnp_gene_grch37_db = work_dir + "reference/dbsnp/dbsnp.gene.grch37.sqlite"

mutate_TIMES_file = work_dir + "mutate_TIMES_py.csv"

cytogenic_refered_file = work_dir + "cytogenic_refered_py.csv"
afterbim_file = work_dir + "afterbim_py.csv"

g1_samples_file = work_dir + "reference/g1_samples.csv"

# plink related define
plink_command_path = 'C:/Users/hykang/Downloads/plink_win64_20210606/plink'
plink_command_option = \
  '--memory 1536 --r2 --ld-window-kb 1000 --ld-window 5000 --ld-window-r2 0.7'

input_file_dir = work_dir + 'reference/Plink'
output_file_dir = work_dir + 'reference/Plink'

plink_ref_path_base = work_dir + 'reference/Plink/ref'

before_checkLDPlink_file = work_dir + "before_checkLDPlink_py.csv"
after_merge_g1_samples_file = work_dir + "after_merge_g1_samples_py.csv"
checkLDPlink_output_file = work_dir + "checkLDPlink_output_py.csv"
after_mutate_beta_zero_file = work_dir + "after_mutate_beta_zero_py.csv"
after_groupby_unfiltered_meta_file = work_dir + "after_groupby_unfiltered_meta_py.csv"
before_mutate_snp_paf_file = work_dir + "before_mutate_snp_paf_py.csv"
unfiltered_meta_file = work_dir + "unfiltered_meta_py.csv"
beta_zero_file = work_dir + "beta_zero_py.csv"
final_meta_file = work_dir + "final_meta_py.csv"
