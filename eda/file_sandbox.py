import pandas as pd
from Bio import SeqIO

def convert_fasta_to_parquet(fasta_file_path):
    # Fasta files do not have a fully strict schema;
    # the below only works with this specific fasta type
    fasta_seq = SeqIO.index(fasta_file_path, "fasta")
    fasta = pd.DataFrame()
    # in this case each sequence is a transcript
    for i, (seq_name, v) in enumerate(fasta_seq.items()):
        if i < 100000:
            print(i)
            dscr_list = v.description.split(' ')
            # dscr_str = v.description.split('description')[1]

            row = {'seq': v.seq, 'id': v.id, 'name': v.name, 'description': v.description,
                   'chromosome': dscr_list[2], 'gene': dscr_list[3], 'gene_biotype': dscr_list[4],
                   'tx_biotype': dscr_list[5],
                   # 'gene_symbol': dscr_list[6], # not always present
                   # 'description': dscr_str # not always present
                   }
            fasta = pd.concat([fasta, pd.DataFrame([row])], ignore_index=True)
        else:
            break
    fasta.to_parquet(fasta_index_path + '_100000_rows.parquet')
    print('done')

def convert_fastq_to_parquet(fastq_file_path):
    fastq = pd.DataFrame()
    for i, rec in enumerate(SeqIO.parse(fastq_file_path, "fastq")):
        if i < 100000:
            row = {'id': rec.id, 'name': rec.name, 'seq': rec.seq,
                   'name': rec.name, 'description': rec.description}
            fastq =  pd.concat([fastq, pd.DataFrame([row])], ignore_index=True)
            print('done processing record' + str(i))
        else:
            break
    fastq['seq'] = fastq['seq'].astype(str)
    fastq.to_parquet(fastq_file_path + '_100000_rows.parquet')
    print('done')

# seq data from experiment, sequence info is split across 4 lines
fastq_file_path = '/media/awmundy/Windows/bio/file_sandbox/Thy_O_CX3pos1_IGO_11991_B_7_S39_L001_R1_001.fastq'
# reference genome, sequence info is split across two lines (description, sequence)
fasta_file_path = '/media/awmundy/Windows/bio/file_sandbox/Homo_sapiens.GRCh38.cdna.all.fa'
# not a text file, produced and used by kallisto to perform psuedualignment
fasta_index_path = '/media/awmundy/Windows/bio/file_sandbox/Homo_sapiens.GRCh38.cdna.all.index'
# contain counts, length, effective legth, and transcripts per million for ~100k transcripts (many are 0)
abundance_path = '/media/awmundy/Windows/bio/ac_thymus/rna_txs/fastq_folders/cx3_neg_young_1/abundance.tsv'

# convert_fastq_to_parquet(fastq_file_path)
# convert_fasta_to_parquet(fasta_file_path)

# fasta = pd.read_parquet(fasta_file_path + '_100000_rows.parquet')
# fastq = pd.read_parquet(fastq_file_path + '_100000_rows.parquet')

abun = pd.read_csv(abundance_path, sep='\t')
assert abun.target_id.duplicated().sum() == 0
set(abun.columns) == set(['target_id', 'length', 'eff_length', 'est_counts', 'tpm'])