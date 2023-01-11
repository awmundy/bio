import subprocess
import os
import shutil
from bio.data_prep.config import cfg, run

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



# can't write directly to external drive with prefetch, need to write to temp drive and then move
home_dir = os.path.expanduser('~') + '/'

temp_download_dir = f'{home_dir}temp_sra_download/'
temp_output_dir = f'{temp_download_dir}{run}/rna_txs/'
temp_log_dir = f'{temp_output_dir}fastq_download_logs/'
output_dir = cfg['rna_txs_dir']
os.makedirs(temp_download_dir, exist_ok=True)
os.makedirs(temp_log_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

sra_ids = cfg['sample_map'].keys()

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
os.rmdir(temp_download_dir)