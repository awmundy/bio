import os
import shutil
from datetime import datetime as dt
import subprocess

def rename_fastq_file(fastq_root_dir, fname):
    suffix = fname.split('.fastq')[-1].split('.gz')[0]
    pre_dot_fastq_fname = fname.split('.fastq')[0]
    new_fname = pre_dot_fastq_fname + suffix + '.fastq.gz'

    if fname != new_fname:
        print('renaming', fname, 'to', new_fname)
        os.rename(fastq_root_dir + fname, fastq_root_dir + new_fname)

    return new_fname

def rename_and_move_fastq_files(fastq_root_dir):
    files = [x for x in os.listdir(fastq_root_dir) if x.endswith('.gz')]
    if len(files) == 0:
        print('All fastq files already moved')
        return

    for fname in files:
        # for fastq files with weird extensions like .fastq-004.gz, move "-004" in to filename
        new_fname = rename_fastq_file(fastq_root_dir, fname)

        # store each fastq in its own subdir
        subdir = fastq_root_dir + new_fname.replace('.fastq.gz', '') + '/'
        print('moving file', new_fname, 'to', subdir)
        os.makedirs(subdir)
        shutil.move(fastq_root_dir + new_fname, subdir + new_fname)


def get_fastq_fpaths(fastq_dir):
    """returns a list of all .fastq.gz files in the given fastq_dir's subdirectories"""
    fastq_paths = []
    for subdir in [x.path for x in os.scandir(fastq_dir) if x.is_dir()]:
        fastq_paths += [x.path for x in os.scandir(subdir) if x.path.endswith('.fastq.gz')]

    return fastq_paths


def fastqc(fastq_dir, threads):
    """use fastqc to construct an html summary file of some quality control checks on raw seq data"""

    fastqc_dir = fastq_dir + '/fastqc/'
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)

    fastq_paths = get_fastq_fpaths(fastq_dir)

    # dont overwrite if they already exist
    for _path in fastq_paths:
        fastqc_path = fastqc_dir + _path.split('/')[-1:][0].replace('.fastq.gz', '') + '_fastqc.html'
        if os.path.exists(fastqc_path):
            fastq_paths = [x for x in fastq_paths if x != _path]


    if len(fastq_paths) > 0:
        print('building fastqc output for the following fastq files:', fastq_paths)
        cmd = 'fastqc ' + ' '.join(fastq_paths)
        cmd += f' -t {threads} -o {fastqc_dir}'
        subprocess.run(cmd, shell=True)
    else:
        print('All fastqc output already produced')


def get_index_fpath_from_fasta_fpath(fasta_fpath):
    index_fpath = fasta_fpath.replace('.fa.gz', '') + '.index'
    return index_fpath


def kallisto_build_index(fasta_fpath):
    # first step is to build an index on the reference fasta file.
    #  this will speed up allignment
    index_fpath = get_index_fpath_from_fasta_fpath(fasta_fpath)

    if os.path.exists(index_fpath):
        print('index already exists for', fasta_fpath)
        return

    subprocess.run(f'kallisto index -i {index_fpath} {fasta_fpath}',  shell=True)

    return index_fpath


def kallisto_quant(index_fpath, fastq_dir, threads):

    fastq_fpaths = get_fastq_fpaths(fastq_dir)[:2]
    # todo remove
    fastq_fpaths = fastq_fpaths[:2]
    output_dir = fastq_dir + 'kallisto_quant' + '/'
    os.makedirs(output_dir, exist_ok=True)

    # todo add to args(probably a new arg that is a dict called kallisto_quant_args)
    # number of reads of the polymer that occurred during gene sequencing,
    #  pass --single for single read or empty string for paired (default execution is paired)
    read_end_type = '--single'
    # sequencing procedure has an average fragment length of the sequences (and standard deviation)
    frag_length = 250
    frag_length_sd = 30

    # todo consider not using a log this way and instead capturing stdout manually so that python console
    #   shows it too
    log_file_path = output_dir + 'log_' + dt.now().strftime('%Y_%m_%d_%Hh%Mm') + '.log'

    cmd = (f'kallisto quant '
           f'-i {index_fpath} '
           f'-o {output_dir} '
           f'-t {threads} '
           f'{read_end_type} '
           f'-l {frag_length} '
           f'-s {frag_length_sd} '
           )
    cmd += + ' '.join(fastq_fpaths)
    cmd += f'&> {log_file_path}'
    print(cmd)
    os.system(cmd)

## Flow ##
# <Download fasta file> -> fasta file -> <kallisto> -> index
# <Download fasta files> -> <Move fastq files> -> fastq files -> <fastqc> -> HTML output

# http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
fasta_fpath = '/media/amundy/Windows/diyt/data/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
# files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
fastq_dir = '/media/amundy/Windows/diyt/data/fastq/'

threads = 10

index_fpath = kallisto_build_index(fasta_fpath)
rename_and_move_fastq_files(fastq_dir)
fastqc(fastq_dir, threads)


