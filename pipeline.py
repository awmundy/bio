import os
import shutil
from datetime import datetime as dt
import subprocess
from subprocess import PIPE, STDOUT
import json
import glob

def fix_fastq_gz_extension_if_needed(fastq_root_dir, fname):
    suffix = fname.split('.fastq')[-1].split('.gz')[0]
    pre_dot_fastq_fname = fname.split('.fastq')[0]
    new_fname = pre_dot_fastq_fname + suffix + '.fastq.gz'

    if fname != new_fname:
        print('renaming', fname, 'to', new_fname)
        os.rename(fastq_root_dir + fname, fastq_root_dir + new_fname)

    return new_fname


def assert_one_fastq_gz_file_in_dir(_dir):
    potential_fastq_files = []
    for _file in os.listdir(_dir):
        if (_file.endswith('.gz')) & ('fastq' in _file):
            potential_fastq_files += [_file]
    assert len(potential_fastq_files) == 1

def rename_fastq_files_and_store_each_in_own_subdir(fastq_root_dir):
    '''Moves each fastq into its own subdir, makes a subdir for storing all of these subdirs, and makes
    the fastq files/dirs have consistent file extensions
    params:
        fastq_root_dir: directory where fastq .gz files are stored together
        fastq_folders_dir: directory where fastq .gz files are to be moved, with each of them getting
                          their own subdir
    '''
    # make a directory for storing fastq folders
    fastq_folders_dir = fastq_root_dir + 'fastq_folders/'
    os.makedirs(fastq_folders_dir, exist_ok=True)

    files = get_fastq_fpaths(fastq_root_dir, ignore_fastq_folders_dir=True)
    if len(files) == 0:
        print('All fastq files already moved to their own subdir in', fastq_folders_dir)
        return

    for fpath in files:
        fdir = os.path.dirname(fpath)+ '/'
        fname = os.path.basename(fpath)
        # for fastq files with weird extensions like .fastq-004.gz, move "-004" in to filename
        fname = fix_fastq_gz_extension_if_needed(fastq_root_dir, fname)
        fpath = fdir + fname

        # store each fastq in its own subdir
        subdir = fastq_folders_dir + fname.replace('.fastq.gz', '') + '/'
        print('moving file', fname, 'to', subdir)
        os.makedirs(subdir)
        shutil.move(fpath, subdir + fname)

    for fastq_dir_name in os.listdir(fastq_folders_dir):
        assert_one_fastq_gz_file_in_dir(fastq_folders_dir + fastq_dir_name)


def get_fastq_fpaths(_dir, ignore_fastq_folders_dir=False):
    """returns a list of all .fastq.gz files in the given _dir's subdirectories"""
    # fastq_paths = []
    # for subdir in [x.path for x in os.scandir(_dir) if x.is_dir()]:
    #     fastq_paths += [x.path for x in os.scandir(subdir) if x.path.endswith('.fastq.gz')]
    fastq_paths = glob.glob(_dir + '/**/*.fastq.gz', recursive=True)
    if ignore_fastq_folders_dir:
        fastq_paths = [x for x in fastq_paths if 'fastq_folders' not in x]

    return fastq_paths


def fastqc(fastq_folders_dir, threads):
    """use fastqc to do some quality control checks on raw seq data
    outputs:
        <fastq_filename>_fastqc.html file: html report showing the qc results
        <fastq_filename>_fastqc.zip file: collection of qc results and other inputs to the html report
        """

    fastq_paths = get_fastq_fpaths(fastq_folders_dir)

    # dont overwrite if they already exist
    fastq_paths_where_fastqc_is_needed = []
    for fastq_path in fastq_paths:
        fastq_dir = os.path.dirname(fastq_path) + '/'
        fastq_file_name = os.path.basename(fastq_path)
        fastqc_html_file_name = fastq_file_name.replace('.fastq.gz', '_fastqc.html')
        # don't overwrite
        if not os.path.exists(fastq_dir + fastqc_html_file_name):
            fastq_paths_where_fastqc_is_needed += [fastq_path]

    if len(fastq_paths_where_fastqc_is_needed) > 0:
        print('building fastqc output for the following fastq files:', fastq_paths_where_fastqc_is_needed)
        cmd = 'fastqc ' + ' '.join(fastq_paths_where_fastqc_is_needed)
        # each thread handles one fastq file and writes the output in the dir of the fastq file being processed
        cmd += f' -t {threads}'
        subprocess.run(cmd, shell=True)
    else:
        print('All fastqc output already produced for fastq files in', fastq_folders_dir)

def get_index_fpath_from_ref_genome_fpath(ref_genome_fpath):
    index_fpath = ref_genome_fpath.replace('.fa.gz', '') + '.index'
    return index_fpath


def kallisto_build_index(ref_genome_fpath):

    index_fpath = get_index_fpath_from_ref_genome_fpath(ref_genome_fpath)

    if os.path.exists(index_fpath):
        print('Index already exists for', ref_genome_fpath)
        return index_fpath

    subprocess.run( f'kallisto index -i {index_fpath} {ref_genome_fpath}', shell=True)

    return index_fpath


def kallisto_quant(index_fpath, fastq_folders_dir, threads, seq_params):
    """
    Runs the kallisto quantification process on each fastq file in fastq_root_dir
    params:
        threads: number of threads to use for EACH fastq file being processed. Multiple fastq files
                 are not multiprocessed; fastq files are iterated serially and multiprocessing is
                 used during the quantification algorithm
        seq_params: dict of params describing some sequencing charcteristics of the fastq files
                    that kallisto quant needs to know
    outputs:
        abundance.h5: a HDF5 binary file containing run info, abundance esimates, bootstrap estimates,
                      and transcript length information length. This file can be read in by sleuth in R
        abundance.tsv: a plaintext file of the abundance estimates. It does not contains bootstrap estimates
        run_info.json: metadata about the kallisto quant run
    """

    fastq_paths = get_fastq_fpaths(fastq_folders_dir)

    # build list of fastq paths where kallisto quant hasnt been run yet
    fastq_paths_where_kallisto_quant_is_needed = []
    for fastq_path in fastq_paths:
        fastq_dir = os.path.dirname(fastq_path) + '/'
        if not os.path.exists(fastq_dir + 'abundance.h5'):
            fastq_paths_where_kallisto_quant_is_needed += [fastq_path]

    if len(fastq_paths_where_kallisto_quant_is_needed) == 0:
        print('kallisto quant already complete for all fastq files in', fastq_folders_dir)
        return

    for fastq_path in fastq_paths_where_kallisto_quant_is_needed:
        print('running kallisto quant on', fastq_path)

        # save the output next to the fastq files
        quant_out_dir = os.path.dirname(fastq_path) + '/'

        cmd = (f"kallisto quant "
               f"-i {index_fpath} "
               f"-o {quant_out_dir} "
               f"-t {threads} "
               f"-l {seq_params['frag_length']} "
               f"-s {seq_params['frag_length_sd']} "
               f"{seq_params['read_end_type']} "
               f"{fastq_path}")

        log_path = quant_out_dir + 'kallisto_quant_log.log'
        # todo figure out how to print stuff out as it's running instead of once subprocess has completed
        with open(log_path, 'wb') as f:
            # capture both the stdout and stderr in one object
            output = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=STDOUT)

            # print and write to log
            for line in output.stdout.splitlines():
                print(line)
                f.write(b'\n')
                f.write(line)

            if output.stderr:
                for line in output.error.splitlines():
                    print(line)
                    f.write(b'\n')
                    f.write(line)
                raise Exception('Kallisto quant process failed')
    print('done running kallisto quant')


def multiqc(fastq_root_dir):
    '''multiqc looks for all relevant files in fastq_dir (i.e. all fastqc files, kallisto quant logs) and
    builds a summary report out of them'''
    multiqc_dir = fastq_root_dir + 'multiqc/'
    os.makedirs(multiqc_dir, exist_ok=True)

    # if there are existing multiqc report, multiqc automatically will make new subdirs
    # to save new runs in
    cmd = f'multiqc -d {fastq_root_dir} -o {multiqc_dir}'
    subprocess.run(cmd, shell=True)



## Flow ##
# <Download fasta file> -> fasta file -> <kallisto> -> index
# <Download fastq files> -> <Move fastq files> -> fastq files -> <fastqc> -> HTML output
# index + fastq file -> kallisto quant
# fastqc output + kallisto quant output -> multiqc

cfgs = \
    {'diyt':
        {# http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
         'ref_genome': '/media/amundy/Windows/bio/eference_genomes/human/Homo_sapiens.GRCh38.cdna.all.fa.gz',
         # files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
         'rna_txs': '/media/amundy/Windows/bio/diyt/rna_txs/',
         'seq_params': {'read_end_type': '--single', 'frag_length': 250, 'frag_length_sd': 30}},
    'ac_thymus':
        {# http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/
         'ref_genome': '/media/amundy/Windows/bio/reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
         'rna_txs': '/media/amundy/Windows/bio/ac_thymus/rna_txs/',
         'seq_params': {'read_end_type': '--single', 'frag_length': 250, 'frag_length_sd': 30}}
        }

run_type = 'ac_thymus'
cfg = cfgs[run_type]
threads = 10

index_fpath = kallisto_build_index(cfg['ref_genome'])
rename_fastq_files_and_store_each_in_own_subdir(cfg['rna_txs'])
# fastqc(fastq_folders_dir, threads)
# kallisto_quant(index_fpath, fastq_folders_dir, threads, seq_params)
# multiqc(fastq_root_dir)

