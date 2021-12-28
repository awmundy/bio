import wget
import os
os.environ['R_HOME'] = '/home/awmundy/anaconda3/envs/bio/lib/R' # to get R imports to work
from IPython.display import Image
import rpy2.robjects as robjects
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from Bio import Entrez, SeqIO, Seq
from Bio.Alphabet import IUPAC
import gzip

warnings.simplefilter(action='ignore', category=FutureWarning)


# console needs linux path root (used to read in stuff in jupyer nb)
bio_drive = '/mnt/c/Drive/Bio/'

# can write to c by pointing to the mount in linux
# wget.download('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/historical_data/former_toplevel/sequence.index',
#               '/mnt/c/drive/bio/data/sequence.index')


# some sort of metadata file where each record is metadata for a sequence
t = pd.read_csv('data/sequence.index', delimiter=r'\t', nrows=10, engine='python')
g = t[['RUN_NAME', 'READ_COUNT' ]].groupby(by='RUN_NAME', as_index=False).agg(sum)

g.plot.bar(x="RUN_NAME", y="READ_COUNT", rot=70, title="Read counts by Run")
plt.show(block=True) #:)



Entrez.email = "awmundy@email.here"
hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='fasta')
# get a sequence object
seq = SeqIO.read(hdl, 'fasta')

rna = seq.translate()