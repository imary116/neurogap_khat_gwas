__author__ = 'Mary T. Yohannes'

# This script reruns LD pruning and produces new set of PLINK files using PLINK - input = filtered and LD pruned NeuroGAP PLINK files

import hailtop.batch as hb
import hail as hl

# get file size to allocate resources for batch jobs accordingly
def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

# Step 2 AGAIN: prune to remove highly correlated SNPs
# perform pruning with a window size of 200 variants, sliding across the genome with step size of 100 variants at a time, and filter out any SNPs with LD r^2 > 0.2
def ld_prune(b, bg_files, storage_size):
    j = b.new_job(name='LD prune AGAIN')  # define job and label it #f'{label}-clumping'
    j.image('hailgenetics/genetics:0.2.37')  # Docker image that contains PLINK
    j.storage(f'{storage_size}Gi')  # increase storage according to file size
    j.cpu(16) # set cpu
    j.memory('highmem')

    j.declare_resource_group(ofile={
        'prune.in': '{root}.prune.in',
        'rmdup.list': '{root}.rmdup.list',
        'log': '{root}.log'
    })

    j.command(f'''
    plink2 --bfile {bg_files} \
    --rm-dup force-first list \
    --indep-pairwise 200 100 0.2 \
    --out {j.ofile} ''')

    j.command(f"echo {j.ofile['rmdup.list']}") # in order to write it out by itself

    return j

# Step 3: generate a new set of PLINK files with filtered variants and samples
def generate_files(b, plink_files, pruned_snps, storage_size):
    j = b.new_job(name='generate a filtered set of PLINK files AGAIN')  # define job and label it
    j.image('hailgenetics/genetics:0.2.37')  # Docker image that contains PLINK
    j.storage(f'{storage_size}Gi')  # increase storage according to file size
    j.cpu(16)
    j.memory('highmem')

    j.declare_resource_group(ofile={
        'bed': '{root}.bed',
        'bim': '{root}.bim',
        'fam': '{root}.fam',
        'log': '{root}.log'})

    j.command(f''' 
    plink2 --bfile {plink_files} \
    --extract {pruned_snps} \
    --make-bed \
    --out {j.ofile} ''')

    return j

if __name__ == '__main__':

    backend = hb.ServiceBackend(billing_project='neale-pumas-bge', remote_tmpdir='gs://neurogap-bge-imputed-regional/mary/khat_gwas/tmp')  # set up backend

    b = hb.Batch(backend=backend, name='perform LD pruning and subset PLINK files AGAIN')  # define batch

    plink_files = b.read_input_group(**{
        'bed': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run1.bed',
        'bim': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run1.bim',
        'fam': 'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run1.fam'})

    # get file size
    plink_size = get_file_size('gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run1.bed')
    storage_size1 = round(plink_size + 10)

    # run Step 2
    run_ld_prune = ld_prune(b, plink_files, storage_size1)
    b.write_output(run_ld_prune.ofile, f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run2')


    # get file size
    storage_size3 = round(plink_size + (plink_size/2) + 5)

    # run Step 3
    run_generate_files = generate_files(b, plink_files, run_ld_prune.ofile['prune.in'], storage_size3)

    # write out filtered and pruned PLINK files
    b.write_output(run_generate_files.ofile, f'gs://neurogap-bge-imputed-regional/mary/khat_gwas/intermediate_files/NeuroGAP_filtered_ldpruned_run2')

    b.run(open=True, wait=False)  # run batch

    backend.close()










