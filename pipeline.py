import os
import shutil
from datetime import datetime as dt
import subprocess
from subprocess import PIPE, STDOUT
import json
import glob

def fix_fastq_file_names_if_needed(fpaths):

    new_fpaths = []

    for fpath in fpaths:
        fdir = os.path.dirname(fpath)+ '/'
        fname = os.path.basename(fpath)
        
        # for fastq files with weird extensions like .fastq-004.gz, move "-004" in to filename
        suffix = fname.split('.fastq')[-1].split('.gz')[0]
        pre_dot_fastq_fname = fname.split('.fastq')[0]
        new_fname = pre_dot_fastq_fname + suffix + '.fastq.gz'

        if fname != new_fname:
            # rename file
            print('renaming', fname, 'to', new_fname)
            new_fpath = fdir + new_fname
            os.rename(fpath, new_fpath)
            new_fpaths += [new_fpath]
        else:
            # or just append current path to list to return
            new_fpaths += [fpath]

    return new_fpaths


def assert_one_fastq_gz_file_in_dir(_dir):
    potential_fastq_files = []
    for _file in os.listdir(_dir):
        if (_file.endswith('.gz')) & ('fastq' in _file):
            potential_fastq_files += [_file]
    assert len(potential_fastq_files) == 1

def move_each_fastq_file_to_unique_dir(fpaths, fastq_folders_dir):
    for fpath in fpaths:
        fname = os.path.basename(fpath)
        subdir = fastq_folders_dir + fname.replace('.fastq.gz', '') + '/'
        print('moving file', fname, 'to', subdir)
        os.makedirs(subdir)
        shutil.move(fpath, subdir + fname)

def prep_fastq_files_diyt(rna_txs_dir, fastq_folders_dir):
    fpaths = get_fastq_fpaths(rna_txs_dir, ignore_fastq_folders_dir=True)
    if len(fpaths) == 0:
        print('All fastq files already moved to their own subdir in', fastq_folders_dir)
        return

    fpaths = fix_fastq_file_names_if_needed(fpaths)
    move_each_fastq_file_to_unique_dir(fpaths, fastq_folders_dir)

    for fastq_dir_name in os.listdir(fastq_folders_dir):
        assert_one_fastq_gz_file_in_dir(fastq_folders_dir + fastq_dir_name)

def prep_fastq_files(cfg, run_type):
    '''Moves each fastq into its own subdir, makes a subdir for storing all of these subdirs, and makes
    the fastq files/dirs have consistent file extensions
    params:
        rna_txs_dir: directory where fastq .gz files are stored together
        fastq_folders_dir: directory where fastq .gz files are to be moved, with each of them getting
                          their own subdir
    '''
    rna_txs_dir = cfg['rna_txs_dir']

    # make a directory for storing fastq folders
    fastq_folders_dir = rna_txs_dir + 'fastq_folders/'
    os.makedirs(fastq_folders_dir, exist_ok=True)

    if run_type == 'ac_thymus':
        raw_rna_txs_dir = cfg['raw_rna_txs_dir']
        prep_fastq_files_ac_thymus(raw_rna_txs_dir, fastq_folders_dir)

    elif run_type == 'diyt':
        prep_fastq_files_diyt(rna_txs_dir, fastq_folders_dir)

    else:
        raise Exception('Invalid run_type')

def get_fastq_fpaths(_dir, ignore_fastq_folders_dir=False):
    """returns a list of all .fastq.gz files in the given _dir's subdirectories"""
    # fastq_paths = []
    # for subdir in [x.path for x in os.scandir(_dir) if x.is_dir()]:
    #     fastq_paths += [x.path for x in os.scandir(subdir) if x.path.endswith('.fastq.gz')]
    fastq_paths = glob.glob(_dir + '/**/*.fastq.gz', recursive=True)
    if ignore_fastq_folders_dir:
        fastq_paths = [x for x in fastq_paths if 'fastq_folders' not in x]

    return fastq_paths


def fastqc(rna_txs_dir, threads):
    """use fastqc to do some quality control checks on raw seq data
    outputs:
        <fastq_filename>_fastqc.html file: html report showing the qc results
        <fastq_filename>_fastqc.zip file: collection of qc results and other inputs to the html report
        """
    fastq_folders_dir = rna_txs_dir + 'fastq_folders/'
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

def get_fastq_files_list_and_assert_n_files(sample_fastq_folder, expected_number_of_fastq_files):
    fastq_files = [sample_fastq_folder + x for x in os.listdir(sample_fastq_folder) if x.endswith('.fastq.gz')]
    assert len(fastq_files) == expected_number_of_fastq_files

    return fastq_files


def kallisto_quant(index_fpath, rna_txs_dir, threads, seq_params):
    """
    Runs the kallisto quantification process on each fastq file (or files for paired reads) in each fastq folder
    in the fastq_folders dir of rna_txs_dir
    params:
        threads: number of threads to use for EACH fastq file being processed. Multiple fastq files
                 are not multiprocessed; fastq files are iterated serially and multiprocessing is
                 used during the quantification algorithm
        seq_params: dict of params describing some sequencing charcteristics of the fastq files
                    that kallisto quant needs to know
    outputs:
        abundance.h5: a HDF5 binary file containing run info, abundance esimates, bootstrap estimates,
                      and transcript length information. This file can be read in by sleuth in R
        abundance.tsv: a plaintext file of the abundance estimates. It does not contains bootstrap estimates
        run_info.json: metadata about the kallisto quant run
    """
    fastq_folders_dir = rna_txs_dir + 'fastq_folders/'

    # build list of fastq paths where kallisto quant hasnt been run yet
    sample_fastq_folders_to_quant = []
    sample_fastq_folder_labels = os.listdir(fastq_folders_dir)
    for _dir in sample_fastq_folder_labels:
        sample_fastq_folder = fastq_folders_dir + _dir + '/'
        if not os.path.exists(sample_fastq_folder + 'abundance.h5'):
            sample_fastq_folders_to_quant += [sample_fastq_folder]


    if len(sample_fastq_folders_to_quant) == 0:
        print('kallisto quant already complete for all fastq folders in', fastq_folders_dir)
        return

    for sample_fastq_folder in sample_fastq_folders_to_quant:
        print('running kallisto quant on', sample_fastq_folder)

        cmd = (f"kallisto quant "
               f"-i {index_fpath} "
               f"-o {sample_fastq_folder} " # output directory == fastq file(s) location
               f"-t {threads} ")

        if seq_params['read_end_type'] == '--single':
            fastq_files = get_fastq_files_list_and_assert_n_files(sample_fastq_folder, 1)
            fastq_file_path = fastq_files[0]
            cmd += seq_params['read_end_type'] + ' '
            cmd += f"-l {seq_params['frag_length']} "
            cmd += f"-s {seq_params['frag_length_sd']} "
            cmd += fastq_file_path
        elif seq_params['read_end_type'] == '--double':
            fastq_files = get_fastq_files_list_and_assert_n_files(sample_fastq_folder, 2)
            cmd += ' '.join(fastq_files)

        log_path = sample_fastq_folder + 'kallisto_quant_log.log'
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

def multiqc(rna_txs_dir):
    '''multiqc looks for all relevant files in fastq_dir (i.e. all fastqc files, kallisto quant logs) and
    builds a summary report out of them'''
    multiqc_dir = rna_txs_dir + 'multiqc/'
    os.makedirs(multiqc_dir, exist_ok=True)

    # if there are existing multiqc report, multiqc automatically will make new subdirs
    # to save new runs in
    cmd = f'multiqc -d {rna_txs_dir} -o {multiqc_dir}'
    subprocess.run(cmd, shell=True)

def get_single_lane_fastq_file_path(raw_sample_dir, lane_label, read_type):
    file_path_candidates = [raw_sample_dir + x for x in os.listdir(raw_sample_dir) if lane_label + read_type in x]
    assert len(file_path_candidates) == 1
    file_path = file_path_candidates[0]
    assert os.path.exists(file_path)
    return file_path

def get_ac_thymus_sample_rename_map():
    sample_rename = {'Sample_Thy_O_CX3pos1_IGO_11991_B_7': 'cx3_pos_old_1',
                     'Sample_Thy_O_CX3pos2_IGO_11991_B_9': 'cx3_pos_old_2',
                     'Sample_Thy_O_CX3pos3_IGO_11991_B_11': 'cx3_pos_old_3',
                     'Sample_Thy_Y_CX3neg1_IGO_11991_B_2': 'cx3_neg_young_1',
                     'Sample_Thy_Y_CX3neg2_IGO_11991_B_4': 'cx3_neg_young_2',
                     'Sample_Thy_Y_CX3neg3_IGO_11991_B_6': 'cx3_neg_young_3',
                     'Sample_Thy_Y_CX3pos1_IGO_11991_B_1': 'cx3_pos_young_1',
                     'Sample_Thy_Y_CX3pos2_IGO_11991_B_3': 'cx3_pos_young_2'
                     }
    return sample_rename

def prep_fastq_files_ac_thymus(raw_rna_txs_dir, fastq_folders_dir):

    sample_rename = get_ac_thymus_sample_rename_map()

    # create dict mapping raw file dir to combined file dir
    dirs_to_process = {}
    for raw_sample_label in os.listdir(raw_rna_txs_dir):
        # get new label for the sample folders and files
        new_sample_label = sample_rename[raw_sample_label]
        combined_sample_dir = fastq_folders_dir + new_sample_label + '/'
        if not os.path.exists(combined_sample_dir):
            dirs_to_process[new_sample_label] = {'raw_sample_dir': raw_rna_txs_dir + raw_sample_label + '/',
                                                 'combined_sample_dir': combined_sample_dir}
    if (len(dirs_to_process) == 0):
        print('no combining of fastq files necessary')
        return

    for new_sample_label, raw_and_combined_dirs in dirs_to_process.items():
        for read_type in ['R1', 'R2']:
            raw_sample_dir = raw_and_combined_dirs['raw_sample_dir']
            combined_sample_dir = raw_and_combined_dirs['combined_sample_dir']

            fpath_1 = get_single_lane_fastq_file_path(raw_sample_dir, 'L001_', read_type)
            fpath_2 = get_single_lane_fastq_file_path(raw_sample_dir, 'L002_', read_type)

            os.makedirs(combined_sample_dir, exist_ok=True)
            combined_fpath = combined_sample_dir + new_sample_label + '_' + read_type +'.fastq.gz'

            cmd = f'cat {fpath_1} {fpath_2} > {combined_fpath}'
            subprocess.run(cmd, shell=True)

            assert os.path.getsize(combined_fpath) > 0
            print(f'done combining files {fpath_1} {fpath_2}')


## Flow ##
# <Download fasta file> -> fasta file -> <kallisto> -> index
# <Download fastq files> -> <Move fastq files> -> fastq files -> <fastqc> - > HTML output
# index + fastq file -> kallisto quant
# fastqc output + kallisto quant output -> multiqc

cfgs = \
    {'diyt':
        {# http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
         'ref_genome': '/media/amundy/Windows/bio/reference_genomes/human/Homo_sapiens.GRCh38.cdna.all.fa.gz',
         # files source (course dataset): https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7
         'rna_txs_dir': '/media/amundy/Windows/bio/diyt/rna_txs/',
         'seq_params': {'read_end_type': '--single', 'frag_length': 250, 'frag_length_sd': 30}},
    'ac_thymus':
        {# http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/
         'ref_genome': '/media/amundy/Windows/bio/reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
         'raw_rna_txs_dir': '/media/amundy/Windows/bio/ac_thymus/raw_rna_txs/' ,
         'rna_txs_dir': '/media/amundy/Windows/bio/ac_thymus/rna_txs/',
          # todo figure out these (and other) params for the thymus data
         'seq_params': {'read_end_type': '--double'}}
        }

run_type = 'ac_thymus'
cfg = cfgs[run_type]
threads = 10

index_fpath = kallisto_build_index(cfg['ref_genome'])
prep_fastq_files(cfg, run_type)
fastqc(cfg['rna_txs_dir'], threads)
kallisto_quant(index_fpath, cfg['rna_txs_dir'], threads, cfg['seq_params'])
multiqc(cfg['rna_txs_dir'])

# Thy_Y_CX3pos3_IGO_11991_B_5_S37_R1_001.fastq.gz
#   - not in gzip format
#   - kallisto quant fails to pseudoalign
# Thy_Y_CX3pos3_IGO_11991_B_5_S37_R2_001.fastq.gz
#   - 45% of way through fastqc: uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Midline
#     <seq string> didn't start with '+'
#   - causes kallisto quant to hang?