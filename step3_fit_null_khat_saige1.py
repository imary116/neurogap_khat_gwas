__author__ = 'Mary T. Yohannes'

# This script fits the null logistic/linear mixed model using a full GRM using SAIGE
## https://saigegit.github.io/SAIGE-doc/docs/single_step1.html

import hail as hl
import hailtop.batch as hb

# get file size to allocate resources for batch jobs accordingly
def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

# fit the null logistic/linear mixed model using a full GRM
def fit_null_model(b, site, plink_files, phenotype_file, phenotype_col, covariate_cols, trait_type, invNormalize, storage_size):
    j = b.new_job(name=f'{site} - {phenotype_col}')
    j.image('wzhou88/saige:1.3.0') # SAIGE docker image
    j.cpu(16)
    j.storage(f'{storage_size}Gi')  # increase storage according to file size

    # model and variance ratio files - inputs for step 4 single-variant association tests (SAIGE step 2 - used in step4_association_test_khat_saige2.py script)
    j.declare_resource_group(ofile={
        'rda': '{root}.rda',
        'varianceRatio.txt': '{root}.varianceRatio.txt'})

    j.command(f'''
    step1_fitNULLGLMM.R \
    --plinkFile={plink_files} \
    --phenoFile={phenotype_file} \
    --phenoCol={phenotype_col} \
    --covarColList={covariate_cols} \
    --sampleIDColinphenoFile=IID \
    --traitType={trait_type} \
    --outputPrefix={j.ofile} \
    --nThreads=16 \
    --LOCO=True \
    --invNormalize={invNormalize} ''')

    return j


if __name__ == '__main__':

    backend = hb.ServiceBackend(billing_project='neale-pumas-bge', remote_tmpdir='gs://neurogap-bge-imputed-regional/mary/khat_gwas/tmp')  # set up backend

    # phenotypes we are interested in
    phenotype_cols = ["assist_khat", "assist_khat_amt"]

    for phenotype in phenotype_cols:

        b = hb.Batch(backend=backend, name=f'step3: SAIGE GWAS step1 - {phenotype}')  # define batch

        # study sites
        sites = ["AAU", "KEMRI", "Moi", "UCT", "Uganda"]

        # PLINK files with filtered and pruned (x3) snps - generated using step2_filter_snps.py, step2_ldprune_rerun.py and step2_ldprune_rerun2.py scripts
        plink_files = b.read_input_group(**{
            'bed': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run3.bed',
            'bim': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run3.bim',
            'fam': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run3.fam'})

        for site in sites:

            # phenotype files by site - generated using step1_generate_phenotype_file_by_site.Rmd script
            phenotype_file = b.read_input(f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/phenotype_files/{site}_khat_pheno.txt')

            # get file size
            plink_size = get_file_size('gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run3.bed')
            storage_size = round(plink_size + 3)

            # covariates in the phenotype file
            covariate_cols = "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex"

            # adjust parameters based on what type of trait the phenotype is
            if phenotype == "assist_khat":
                trait_type = "binary"
                invNormalize = "FALSE"
            else:
                trait_type = "quantitative"
                invNormalize = "TRUE"

            # run function
            run_fit_null = fit_null_model(b, site, plink_files, phenotype_file, phenotype, covariate_cols, trait_type, invNormalize, storage_size)

            b.write_output(run_fit_null.ofile, f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/saige_step1_outputs/{phenotype}/{site}_{phenotype}_saige_step1')

        b.run(wait=False)  # run batch

    backend.close()  # close a batch backend
