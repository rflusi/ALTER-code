#####################
# import statements #
#####################

import os
import pandas as pd
import sys
import pysam
import pickle
import gzip
import tqdm
from tabulate import tabulate
from Bio import Seq
import re
from pandarallel import pandarallel
import gc

#############
# functions #
#############

###############
# main script #
###############
if __name__ == '__main__':
    #---Step 1: Define Variables ---
    # from args
    proj_dir = sys.argv[1]
    var_dir = sys.argv[2]
    ref_dir = sys.argv[3]

    # directories
    align_dir = os.path.join(var_dir, 'transcript-align')

    # ref files
    ref_genome_name = os.path.split(ref_dir)[1]
    transcript_fa_path = os.path.join(ref_dir, f'{ref_genome_name}.transcripts.filtered.fa.gz')

    #---Step 2: build dict of transcript contigs ---
    print(f'\nBuilding Transcript Contig Lookup for {ref_genome_name} Transcripts\n')
    
    transcript_contig_path = os.path.join(align_dir, 'transcript-id-to-contig.pkl.gz')
    if os.path.isfile(transcript_contig_path):
        with gzip.open(transcript_contig_path, 'rb') as f:
            transcript_contig_dict = pickle.load(f)
        print(f'\nLoaded Previously Built Transcript Contig Lookup for {len(transcript_contig_dict)} Transcripts.\n')
    else:    
        transcript_contig_dict = {}
        for transcript_contig in pysam.FastaFile(transcript_fa_path).references:
            transcript_id = transcript_contig.strip().split('|')[0]
            transcript_contig_dict[transcript_id] = transcript_contig
        
        with gzip.open(transcript_contig_path, 'wb') as f:
            pickle.dump(transcript_contig_dict, f)
        print(f'\nBuilt Transcript Contig Lookup for {len(transcript_contig_dict)} Transcripts.\n')
