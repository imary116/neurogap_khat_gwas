__author__ = 'Mary T. Yohannes'

# This script performs meta-analysis using PLINK after combining each chromosome file by site for each phenotype
# input files are the outputs from step4_association_test_khat_saige2.py script - 5 sites * 22 chromosomes * 2 phenotypes = 220 files in total

import hail as hl
import hailtop.batch as hb

# run meta-analysis on SAIGE association test outputs
def meta_analysis(b, input_files, phenotype, trait_arg):
    j = b.new_job(name=f'{phenotype} - meta-analysis')
    j.image('hailgenetics/genetics:0.2.37')  # Docker image that contains PLINK
    j.cpu(4)
    j.storage('40Gi')

    # setting up files to be concatenated together by site
    j.command('mkdir io/tmp')  # make temporary directory

    for path in input_files:
        file = b.read_input(path['path']) # read in file
        j.command(f'mv {file} io/tmp')  # move file into the temporary directory

    #j.command('ls io/tmp')  # list files in temporary directory

    # concatenate files by site
    j.command(f''' for site in AAU KEMRI Moi UCT Uganda; do cat io/tmp/"$site"_chr*_{phenotype}_saige_step2 > io/tmp/"$site"_allchr_{phenotype}_saige_step2; done ''')

    j.command('echo plink started running')

    # run meta-analysis
    j.command(f''' 
    plink \
    --meta-analysis io/tmp/AAU_allchr_{phenotype}_saige_step2 io/tmp/Moi_allchr_{phenotype}_saige_step2 io/tmp/KEMRI_allchr_{phenotype}_saige_step2 io/tmp/UCT_allchr_{phenotype}_saige_step2  io/tmp/Uganda_allchr_{phenotype}_saige_step2 + {trait_arg}\
    --meta-analysis-a1-field Allele1 \
    --meta-analysis-a2-field Allele2 \
    --meta-analysis-bp-field POS \
    --meta-analysis-chr-field CHR \
    --meta-analysis-snp-field MarkerID \
    --out io/tmp/saige_all_sites_{phenotype} ''')

    # the plink output file has spacing that's inconsistent with R so fix that here
    j.command(f''' awk -v OFS='\t' '$1=$1' io/tmp/saige_all_sites_{phenotype}.meta > io/tmp/saige_{phenotype}_LOCO.meta ''')

    j.command(f'mv io/tmp/saige_{phenotype}_LOCO.meta {j.output}')

    b.write_output(j.output, f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/meta_analysis_outputs/saige_{phenotype}_LOCO.meta')


if __name__ == '__main__':

    backend = hb.ServiceBackend(billing_project='neale-pumas-bge', remote_tmpdir='gs://neurogap-bge-imputed-regional/mary/khat_gwas/tmp')  # set up backend

    b = hb.Batch(backend=backend, name=f'step5: khat GWAS meta-analysis')  # define batch

    # phenotypes we are interested in
    phenotype_cols = ["assist_khat", "assist_khat_amt"]

    for phenotype in phenotype_cols:

        # get file paths
        input_files = hl.hadoop_ls(f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/saige_step2_outputs/{phenotype}/*_chr*_{phenotype}_saige_step2')

        # adjust arguments based on what type of trait the phenotype is - logscale is for binary and qt is quantitative
        if phenotype == "assist_khat":
            trait_arg = "logscale"
        else:
            trait_arg = "qt"

        # run function
        run_meta_analysis = meta_analysis(b, input_files, phenotype, trait_arg)

        b.run(open=True, wait=False)

    backend.close()
