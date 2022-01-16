import os
import shutil
from datetime import datetime as dt
import subprocess
from subprocess import PIPE, STDOUT
import json

def rename_fastq_file(fastq_root_dir, fname):
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

def rename_fastq_files_and_store_each_in_own_subdir(fastq_root_dir, fastq_folders_dir):
    '''Moves each fastq into its own subdir, makes a subdir for storing all of these subdirs, and makes
    the fastq files/dirs have consistent file extensions
    params:
        fastq_root_dir: directory where fastq .gz files are stored together
        fastq_folders_dir: directory where fastq .gz files are to be moved, with each of them getting
                          their own subdir
    '''
    # make a directory for storing fastq folders
    os.makedirs(fastq_folders_dir, exist_ok=True)

    files = [x for x in os.listdir(fastq_root_dir) if x.endswith('.gz')]
    if len(files) == 0:
        print('All fastq files already moved to their own subdir in', fastq_folders_dir)
        return

    for fname in files:
        # for fastq files with weird extensions like .fastq-004.gz, move "-004" in to filename
        new_fname = rename_fastq_file(fastq_root_dir, fname)

        # store each fastq in its own subdir
        subdir = fastq_root_dir + new_fname.replace('.fastq.gz', '') + '/'
        print('moving file', new_fname, 'to', subdir)
        os.makedirs(subdir)
        shutil.move(fastq_root_dir + new_fname, subdir + new_fname)

    for fastq_dir_name in os.listdir(fastq_folders_dir):
        assert_one_fastq_gz_file_in_dir(fastq_folders_dir + fastq_dir_name)


def get_fastq_fpaths(_dir):
    """returns a list of all .fastq.gz files in the given fastq_dir's subdirectories"""
    fastq_paths = []
    for subdir in [x.path for x in os.scandir(_dir) if x.is_dir()]:
        fastq_paths += [x.path for x in os.scandir(subdir) if x.path.endswith('.fastq.gz')]

    return fastq_paths


def fastqc(fastq_folders_dir, threads):
    """use fastqc to construct an html summary file of some quality control checks on raw seq data"""

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
        print('All fastqc output already produced')

def get_index_fpath_from_fasta_fpath(fasta_fpath):
    index_fpath = fasta_fpath.replace('.fa.gz', '') + '.index'
    return index_fpath


def kallisto_build_index(fasta_fpath):

    index_fpath = get_index_fpath_from_fasta_fpath(fasta_fpath)

    if os.path.exists(index_fpath):
        print('Index already exists for', fasta_fpath)
        return index_fpath

    subprocess.run( f'kallisto index -i {index_fpath} {fasta_fpath}', shell=True)

    return index_fpath



def get_new_fastq_inputs_by_comparing_to_kallisto_json(output_dir, fastq_fpaths):
    # check if there was a previous kallisto quant run done by looking for a run_info.json
    json_meta = output_dir + 'run_info.json'
    if os.path.exists(json_meta):
        with open(json_meta) as f:
            meta = json.load(f)
        previous_fastq_inputs = [x for x in meta['call'].split(' ') if '.fastq.gz' in x]
        new_fastq_inputs = [x for x in fastq_fpaths if x not in previous_fastq_inputs]
    else:
        new_fastq_inputs = fastq_fpaths

    return new_fastq_inputs

def kallisto_quant(index_fpath, fastq_dir, threads):
    # todo dont run if there is already stuff in the output dir

    # todo run each fastq serially, when running them all together it stacks them and just produces one
    #   combined set of outputs that aren't meaningful

    fastq_fpaths = get_fastq_fpaths(fastq_dir)[:1]
    output_dir = fastq_dir + 'kallisto_quant' + '/'
    os.makedirs(output_dir, exist_ok=True)

    new_fastq_inputs = get_new_fastq_inputs_by_comparing_to_kallisto_json(output_dir, fastq_fpaths)
    if len(new_fastq_inputs) == 0:
        print('kallisto quant already run on all fastq files in', fastq_dir)
        return

    # todo add to args(probably a new arg that is a dict called kallisto_quant_args)
    # number of reads of the polymer that occurred during gene sequencing,
    #  pass --single for single read or empty string for paired (default execution is paired)
    read_end_type = '--single'
    # sequencing procedure has an average fragment length of the sequences (and standard deviation)
    frag_length = 250
    frag_length_sd = 30

    cmd = (f'kallisto quant '
           f'-i {index_fpath} '
           f'-o {output_dir} '
           f'-t {threads} '
           f'{read_end_type} '
           f'-l {frag_length} '
           f'-s {frag_length_sd} '
           )
    cmd += ' '.join(fastq_fpaths)

    log_fpath = output_dir + 'log.log'

    # todo figure out how to print stuff out as it's running
    with open(log_fpath, 'wb') as f:
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

    print('done')


def multiqc(fastq_root_dir):
    '''multiqc looks for all relevant files in fastq_dir (i.e. all fastqc files, kallisto quant logs'''
    multiqc_dir = fastq_root_dir + 'multiqc/'
    os.makedirs(multiqc_dir, exist_ok=True)

    cmd = f'multiqc -d {fastq_root_dir} -o {multiqc_dir}'
    subprocess.run(cmd, shell=True)



## Flow ##
# <Download fasta file> -> fasta file -> <kallisto> -> index
# <Download fastq files> -> <Move fastq files> -> fastq files -> <fastqc> -> HTML output

# http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
fasta_fpath = '/media/amundy/Windows/diyt/data/fasta/Homo_sapiens.GRCh38.cdna.all.fa.gz'
# files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
fastq_root_dir = '/media/amundy/Windows/diyt/data/fastq/'
fastq_folders_dir = fastq_root_dir + 'fastq_folders/'
threads = 10


index_fpath = kallisto_build_index(fasta_fpath)
rename_fastq_files_and_store_each_in_own_subdir(fastq_root_dir,fastq_folders_dir)
fastqc(fastq_folders_dir, threads)
# kallisto_quant(index_fpath, fastq_root_dir, threads)
# multiqc(fastq_root_dir)


# todo add hdf5 to conda env




# logging.basicConfig(filename=log_path, level=logging.DEBUG)
# cmd = 'conda env list'
# output = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
# logging.info(output.stdout)

