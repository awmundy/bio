# Location of reference genomes, study design, and downstream R outputs
project_dir = '/media/awmundy/Windows/bio/'
# Flexibility offered to store RNA txs on separate drive for space reasons
rna_txs_root_dir = '/media/awmundy/TOSHIBA EXT/'

cfgs = {'senescence_skeletal_muscle_myofiber':
         {'sample_map':
           {'SRR15931118': 'young_1',
            'SRR15931119': 'young_2',
            'SRR15931120': 'young_3',
            'SRR15931121': 'young_4',
            'SRR15931122': 'old_1',
            'SRR15931123': 'old_2',
            'SRR15931124': 'old_3',
            'SRR15931125': 'old_4',
            'SRR15931126': 'old_senolytics_1',
            'SRR15931127': 'old_senolytics_2',
            'SRR15931128': 'old_senolytics_3',
            'SRR15931129': 'old_senolytics_4',
            },
          'paper': 'https://www.nature.com/articles/s43587-022-00250-8',
          'geo': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172254',
          'ref_genome': f'{project_dir}reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
          'rna_txs_dir': f'{rna_txs_root_dir}senescence_skeletal_muscle_myofiber/rna_txs/',
          'seq_params': {'read_end_type': '--double'}},

        'age_related_steatohepatitis':
         {'sample_map':
           {'SRR22048075': 'old_1',
            'SRR22048076': 'old_2',
            'SRR22048077': 'old_3',
            'SRR22048078': 'young_1',
            'SRR22048079': 'young_2',
            'SRR22048080': 'young_3',},
          'paper': 'https://www.nature.com/articles/s43587-022-00348-z',
          'geo': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216592',
          'ref_genome': f'{project_dir}reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
          'rna_txs_dir': f'{rna_txs_root_dir}age_related_steatohepatitis/rna_txs/',
          'seq_params': {'read_end_type': '--double'}
          },

        'diyt':
         {'sample_map':
           {},
          'rna_txs_dir': f'{rna_txs_root_dir}diyt/rna_txs/',
          # http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/
          'ref_genome': f'{project_dir}reference_genomes/human/Homo_sapiens.GRCh38.cdna.all.fa.gz',
          'seq_params': {'read_end_type': '--single', 'frag_length': 250, 'frag_length_sd': 30}},

        'ac_thymus':
         {'sample_map':
           {'Sample_Thy_O_CX3pos1_IGO_11991_B_7': 'cx3_pos_old_1',
            'Sample_Thy_O_CX3pos2_IGO_11991_B_9': 'cx3_pos_old_2',
            'Sample_Thy_O_CX3pos3_IGO_11991_B_11': 'cx3_pos_old_3',
            'Sample_Thy_Y_CX3neg1_IGO_11991_B_2': 'cx3_neg_young_1',
            'Sample_Thy_Y_CX3neg2_IGO_11991_B_4': 'cx3_neg_young_2',
            'Sample_Thy_Y_CX3neg3_IGO_11991_B_6': 'cx3_neg_young_3',
            'Sample_Thy_Y_CX3pos1_IGO_11991_B_1': 'cx3_pos_young_1',
            'Sample_Thy_Y_CX3pos2_IGO_11991_B_3': 'cx3_pos_young_2'},
          # http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/
          'ref_genome': f'{project_dir}reference_genomes/mouse/Mus_musculus.GRCm39.cdna.all.fa.gz',
          'raw_rna_txs_dir': f'{rna_txs_root_dir}ac_thymus/raw_rna_txs/',
          'rna_txs_dir': f'{rna_txs_root_dir}ac_thymus/rna_txs/',
          'seq_params': {'read_end_type': '--double'}}}

run = 'age_related_steatohepatitis'
cfg = cfgs[run]