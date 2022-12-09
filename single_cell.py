import subprocess
import os

def build_index_and_transcript_to_gene_mapping(cfg):
    '''
    Downloads a pre-built index and builts a transcript to gene mapping. If the index is
    already exists, it should not overwrite. This index
    '''
    # index will be downloaded to this location
    index_path = cfg['index_path']
    # tx to gene mapping will be built and stored at this location
    tx_to_gene_path = cfg['tx_to_gene_path']

    if os.path.exists(index_path) & os.path.exists(tx_to_gene_path):
        print('Genome index and trascript-to-gene mapping already produced, not recreating them')
    else:
        cmd = f"kb ref -i  {cfg['index_path']} \
               -d human \
               -g {cfg['tx_to_gene_path']}"
        subprocess.run(cmd, shell=True)

def build_single_seq_counts(cfg):
    '''
    Build counts for single cell data using kb count. Under the hood,
    kb count uses kallisto to pseudoalign and kb bustools to build counts
    '''
    for _file in [cfg['index_path'], cfg['tx_to_gene_path'],
                  cfg['forward_fastq_path'], cfg['reverse_fastq_path']]:
        if not os.path.exists(_file):
            raise Exception(f"{_file} does not exist")
    if os.path.exists(cfg['kallisto_bus_output_dir']):
        print(f"Counts already exist from previous kallisto bus run, not recreating. \n "
              f"To recreate them, delete {cfg['kallisto_bus_output_dir']}")
    else:
        # construct counts in cellranger compatible format
        cmd = f"kb count \
                -i {cfg['index_path']} \
                -g {cfg['tx_to_gene_path']} \
                -x {cfg['sequencing_technology']} \
                -t {cfg['threads']} \
                -o {cfg['kallisto_bus_output_dir']} \
                --cellranger \
                --overwrite \
                {cfg['forward_fastq_path']} {cfg['reverse_fastq_path']}"
        subprocess.run(cmd, shell=True)

# paths must be in linux filesystem, kb_python can't handle referencing files in the windows one
cfgs = \
    {'diyt':
         {'index_path':
              '/home/awmundy/Documents/single_cell_starting_folder/Homo_sapiens.GRCh38.cdna.all.index',
          'tx_to_gene_path':
              '/home/awmundy/Documents/single_cell_starting_folder/homo_sapien_transcript_to_gene.txt',
          'forward_fastq_path':
              '/home/awmundy/Documents/single_cell_starting_folder/pbmc_1k_v3_S1_mergedLanes_R1.fastq.gz',
          'reverse_fastq_path':
              '/home/awmundy/Documents/single_cell_starting_folder/pbmc_1k_v3_S1_mergedLanes_R2.fastq-002.gz',
          'sequencing_technology':
              '10XV3',
          'threads':
              '2',
          'kallisto_bus_output_dir': '/home/awmundy/Documents/single_cell_starting_folder/kallisto_bus_outputs/'
          },

        'ac_thymus':
            {
             }
        }

cfg = cfgs['diyt']

build_index_and_transcript_to_gene_mapping(cfg)
build_single_seq_counts(cfg)
print('Done')

