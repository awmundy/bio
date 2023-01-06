import subprocess
import os

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



def get_senescence_skeletal_muscle_myofiber_sra_ids():
    # paper: https://www.nature.com/articles/s43587-022-00250-8

    #geo: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172254

    sra_ids = [
        'SRR15931118',
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
        'SRR15931129',
        ]

    return sra_ids

output_dir = '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/rna_txs/'
log_dir = '/media/awmundy/Windows/bio/senescence_skeletal_muscle_myofiber/logs/'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
sra_ids = get_senescence_skeletal_muscle_myofiber_sra_ids()

# requires sra-toolkit command line tool (apt install sra-toolkit)
for sra_id in sra_ids:
    with open(log_dir + 'sra_download_and_fastq_conversion.txt', 'wb') as log:
        print (f'downloading {sra_id}')
        cmd = f'prefetch {sra_id} --output-directory {output_dir}'
        sra_download_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        write_to_log(sra_download_output, log)

        print(f'converting {sra_id} to fastq')
        cmd = f'fastq-dump --outdir {output_dir} --gzip  --readids --read-filter pass --dumpbase ' \
              f'--split-3 --clip {output_dir}{sra_id}/{sra_id}.sra'
        fastq_conversion_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        write_to_log(fastq_conversion_output, log)

        print(f'done downloading and converting {sra_id}')