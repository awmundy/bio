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

    # todo figure out how to run shell commands with correct (wsl) env
    #   - os and subproccess only see the windows conda envs (and cant figure out how to switch envs anyway)
    #   - try something from here? https://stackoverflow.com/questions/57693460/using-wsl-bash-from-within-python

    gz_files = [x for x in os.listdir(fastq_dir) if x.endswith('.gz')]
    assert len(gz_files) == 1
    fastq_name = gz_files[0]
    # fastq_name = 'SRR8668756.fastq-005'

    fastq_filepath = fastq_dir + fastq_name

    os.system(f'fastqc -t {threads} {fastq_filepath}') # doesnt work bc this points to windows
    cmd = 'wsl.exe -ls'
    os.system(cmd) # works but dont know how to change envs


def kallisto_build_index(fasta_filepath):
    # first step is to build an index on the reference fasta file.
    #  this will speed up allignment
    index_filepath = fasta_filepath + '.index'
    if index_filepath in os.listdir('/'.join(f.split('/')[:-1])):
        print('index already exists for', fasta_filepath)
        return

    os.system(f"kallisto index -i {index_filepath} {fasta_filepath}")


def subprocess_debugging():
    # DID NOT FIGURE OUT HOW TO ACTIVATE ENV ON WSL USING SUBPROCESS
    # exit codes
    # - 1: failed for any reason
    # - 127: command not found
    subprocess.check_call('dir', shell=True)  # pass
    subprocess.check_call('ls', shell=True)  # pass
    subprocess.check_call('wsl', shell=True)  # 127
    subprocess.check_call('wsl.exe', shell=True)  # hangs?
    subprocess.check_call(['wsl.exe', 'conda.exe', 'list'])  # pass but points to windows conda
    subprocess.check_call(['wsl.exe', 'which', 'conda.exe'])  # pass, points to windows conda
    subprocess.check_call(['wsl.exe', 'echo', '$PATH'])  # pass, only shows windows paths, doesnt include linux conda
    subprocess.check_call(['wsl.exe', '/root/anaconda3/bin/conda', 'env', 'list'])  # pass, shows linux envs
    subprocess.check_call('wsl.exe which python', shell=True)  # 1
    subprocess.check_call('wsl.exe where python', shell=True)  # 127
    subprocess.check_call('wsl which python', shell=True)  # 127
    subprocess.check_call('wsl.exe source activate centrifuge', shell=True)  # 2 (can't find /etc/profile/d/conda.sh)

    # bash
    subprocess.Popen(['/bin/bash', '-c', 'a="Apples and oranges" && echo "${a/oranges/grapes}"'])  # pass
    subprocess.Popen(['/bin/bash', '-c', 'a="Apples and oranges" && echo "${a/oranges/grapes}"'])

## Current Setup ##
# - ignoring wsl, just using windows envs
# - this means i can't access wsl directories/executables either i think
# - todo consider writing bash scripts to run fastqc on wsl, call the bash scripts from windows
#       - would need to find a way to access wsl programs like fastqc
# todo consider virtual box with ubuntu and shared drive

## Flow ##
# <Move fastq files> -> fastq files -> <fastqc> -> HTML output
# <Download fasta files> -> fasta files -> <kallisto> -> index

# files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
fastq_root_dir = 'c:/diyt/data/fastq/'
threads = 4
make_fastq_subdirs_and_move_fastq_files(fastq_root_dir)


# Eventually will iterate over list of subdirs dirs
fastq_label = 'SRR8668755.fastq-007'
fastq_subdir = fastq_root_dir + fastq_label + '/'
fastq_filepath = fastq_subdir + fastq_label + '.gz'

## Fastqc ##
# todo not currently working
# fastqc(fastq_root_dir, threads)

# Example command line fastqc call, produces html file
# f'fastqc -t {threads} {fastq_filepath}'


### Kallisto ###
# kallisto is used to allign (pseudoallign technically) a seq to a reference genome
#   in order to see what each section of the seq does (ex: identifies genes in the seq)
# file downloaded from http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/


fasta_dir = 'c:/diyt/data/fasta/'
fasta_file = r'Homo_sapiens.GRCh38.cdna.all.fa'
fasta_filepath = fasta_dir + fasta_file
kallisto_build_index(fasta_filepath)

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

