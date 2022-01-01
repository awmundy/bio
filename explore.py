import os
import shutil
from datetime import datetime as dt
import subprocess

def make_fastq_subdirs_and_move_fastq_files(fastq_root_dir):
    files = [x for x in os.listdir(fastq_root_dir) if x.endswith('.gz')]
    if len(files) == 0:
        print('All fastq files already moved')
        return

    for file in files:
        subdir = fastq_root_dir + file.replace('.gz', '') + '/'
        print('moving file', file, 'to', subdir)
        os.makedirs(subdir)
        shutil.move(fastq_root_dir + file, subdir + file)


def fastqc(fastq_dir, threads):
    """use fastqc to construct an html summary file of some quality control checks on raw seq data"""

    fastqc_dir = fastq_dir + '/fastqc/'
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)
    gz_paths = []
    for subdir in [x.path for x in os.scandir(fastq_dir) if x.is_dir()]:
        gz_paths += [x.path for x in os.scandir(subdir) if x.path.endswith('.gz')]


    cmd = 'fastqc ' + ' '.join(gz_paths)
    cmd += f' -t {threads} -o {fastqc_dir}'
    subprocess.run(cmd, shell=True)


def kallisto_build_index(fasta_filepath):
    # first step is to build an index on the reference fasta file.
    #  this will speed up allignment
    index_filepath = fasta_filepath.replace('.fa.gz', '') + '.index'


    if os.path.exists(index_filepath):
        print('index already exists for', fasta_filepath)
        return


    # fasta_dir = os.path.dirname(fasta_filepath)
    # fasta_dir = '/media/amundy/Windows/diyt/data/fasta/test/'
    # fasta_file = 'Danio_rerio.GRCz11.cdna.all.fa.gz'
    # fasta_filepath = fasta_dir + fasta_file
    #
    # index_file = fasta_file.replace('.fa.gz', '') + '.index'
    # index_filepath = fasta_dir + index_file

    subprocess.run(f'kallisto index -i {index_filepath} {fasta_filepath}',  shell=True)


def kallisto_quant():
    # then we use the quant feature
    # save in subdir labelled with the fastq label
    output_dir = fastq_dir + fastq_label + '/'
    # number of reads of the polymer that occurred during gene sequencing,
    #  pass --single for single read or empty string for paired (default execution is paired)
    read_end_type = '--single'
    # sequencing procedure has an average fragment length of the sequences (and standard deviation)
    frag_length = 250
    frag_length_sd = 30
    log_file_path = fastq_filepath + '_' + dt.now().strftime('%Y_%m_%d_%Hh%Mm') + '.log'
    quant_string = (f'kallisto quant '
                    f'-i {index_filepath} '
                    f'-o {output_dir} '
                    f'-t {threads} '
                    f'{read_end_type} '
                    f'-l {frag_length} '
                    f'-s {frag_length_sd} '
                    f'{fastq_filepath} '
                    f'> {log_file_path}'
                    )
    print(quant_string)
    os.system(quant_string)

## Flow ##
# <Download fasta file> -> fasta file -> <kallisto> -> index
# <Download fasta files> -> <Move fastq files> -> fastq files -> <fastqc> -> HTML output

# http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
fasta_filepath = '/media/amundy/Windows/diyt/data/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
kallisto_build_index(fasta_filepath)

# files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
fastq_dir = '/media/amundy/Windows/diyt/data/fastq/'
make_fastq_subdirs_and_move_fastq_files(fastq_dir)

fastqc(fastq_dir, threads=6)



