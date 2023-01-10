import subprocess
import os
import shutil

def write_to_log(output, log):
    for line in output.stdout.splitlines():
        print(line)
        log.write(b'\n')
        log.write(line)

    if output.stderr:
        for line in output.error.splitlines():
            print(line)
            log.write(b'\n')
            log.write(line)

cfgs = {'senescence_skeletal_muscle_myofiber':
            {'sra_ids':
                 ['SRR15931118',
                  'SRR15931119',
                  'SRR15931120',
                  'SRR15931121',
                  'SRR15931122',
                  'SRR15931123',
                  'SRR15931124',
                  'SRR15931125',
                  'SRR15931126',
                  'SRR15931127',
                  'SRR15931128',
                  'SRR15931129',],
             'paper': 'https://www.nature.com/articles/s43587-022-00250-8',
             'geo': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172254'},
        'age_related_steatohepatitis':
            {'sra_ids':
                ['SRR22048075',
                 'SRR22048076',
                 'SRR22048077',
                 'SRR22048078',
                 'SRR22048079',
                 'SRR22048080',],
             'paper': 'https://www.nature.com/articles/s43587-022-00348-z',
             'geo': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216592'}}

run = 'age_related_steatohepatitis'

# can't write directly to external drive with prefetch, need to write to temp drive and then move
temp_output_dir = f'~/temp_sra_download/{run}/'
temp_rna_txs_dir = f'{temp_output_dir}/rna_txs/'
temp_log_dir = f'{temp_output_dir}/fastq_download_logs/'
output_dir = f'/media/awmundy/TOSHIBA EXT/{run}/'
os.makedirs(temp_rna_txs_dir, exist_ok=True)
os.makedirs(temp_log_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

sra_ids = cfgs[run]['sra_ids']

# requires sra-toolkit command line tool (apt install sra-toolkit)
for sra_id in sra_ids:
    with open(temp_log_dir + 'sra_download_and_fastq_conversion.txt', 'wb') as log:
        print (f'downloading {sra_id}')
        cmd = f'prefetch {sra_id} --output-directory {temp_output_dir}'
        sra_download_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        write_to_log(sra_download_output, log)

        print(f'converting {sra_id} to fastq')
        cmd = f'fastq-dump --outdir {temp_output_dir} --gzip  --readids --read-filter pass --dumpbase ' \
              f'--split-3 --clip {temp_output_dir}{sra_id}/{sra_id}.sra'
        fastq_conversion_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        write_to_log(fastq_conversion_output, log)

        print(f'done downloading and converting {sra_id}')

# move files to storage drive and clean up
shutil.move(temp_output_dir, output_dir)
# os.rmdir(temp_output_dir)