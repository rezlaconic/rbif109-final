import argparse
import re
import requests
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt



'''Script for RBIF-109 Assignment 4 : Homolog Finder

        Script will accept a peptide sequence either in the form of a string or FASTA file and conduct a BLASTP search to 
        gather multiple sequence homologs. Homologs will be binned by percent similarity and 96 homologs across all bins will 
        be chosen and exported to a ready-to-order csv.

        Expected Usage:

        python3 homolog_finder.py --fasta tdt.fasta

        Expected Inputs ( Only 1 required):
        --fasta : path to a fasta file containing amino acid sequence
        --aa_sequence  : string of amino acid sequence 

        Expected Outputs:

        1) 96-well-homologs.csv - CSV assortment of 96-homologs 

        2) blast_alignments.txt - TXT file  containing the formatted alignments from the BLASTP query

        3) blast_results.csv - CSV file containing the analyzed BLASTP records

        4) blast_summary.png - PNG Plots describing the spread of percent similarity of homologs obtained
'''

def generate_blastp_df(sequence,eval_thresh=0.00000000001 ):
    '''Run BLAST on query and generate dataframe'''
    columns = ['Title','Hit ID','Hit Def','Length','Score','Bits','Expect','Num Alignments','Identities','Positives','Gaps','Align Length','Strand','Frame','Query','Query Start','Query End','Match','Subject','Subject Start','Subject End']
    rows = []
    
    #generate blast query
    print(f'Running BLASTP for query {sequence} ...')
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=1000)
    blast_records = NCBIXML.parse(result_handle)

    count = 0
    with open('blast_alignments.txt','w') as file:
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < eval_thresh:
                        count += 1
                        file.write("****Alignment****\n")
                        file.write(f"sequence: {alignment.title}\n")
                        file.write(f"length: {alignment.length} \n")
                        file.write(f"{hsp.query}\n")
                        file.write(f"{hsp.match}\n")
                        file.write(f"{hsp.sbjct}\n")

                        #append to rows
                        rows.append([alignment.title,alignment.hit_id,alignment.hit_def,alignment.length,hsp.score,hsp.bits,hsp.expect,hsp.num_alignments,hsp.identities,hsp.positives,hsp.gaps,hsp.align_length,hsp.strand,hsp.frame,hsp.query,hsp.query_start,hsp.query_end,hsp.match,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end])
                
    print(f"There are {count} similar sequences in Blast output")
    blast_df = pd.DataFrame(rows,columns=columns)

    return blast_df

def main():
    ####### ARG HANDLING #######
    parser = argparse.ArgumentParser()
    parser.add_argument("--aa_sequence", help="aa sequence of protein", type=str)
    parser.add_argument("--fasta",help="FASTA file of aa sequence of interest",type=str) #will read the last FASTA ENTRY

    args = parser.parse_args()

    fasta = args.fasta
    aa_sequence = args.aa_sequence
    seq_name = None

    if fasta:
        print(f'Reading {fasta} ...')
        for seq_record in SeqIO.parse(fasta,"fasta"):
            aa_sequence = seq_record.seq
        print(f'Obtained {aa_sequence} from {fasta}...')

        seq_name = fasta.split('.')[0]
        
    ####### BLASTP SEARCH. EXPORT AS CSV.

    blast_df = generate_blastp_df(aa_sequence)
    blast_df['%Similarity'] = [ float(x/(len(aa_sequence)))*100 for x in list(blast_df['Identities'])]
    blast_df['Scientific Name'] = [re.findall('\[.*?\]',header)[-1].replace('[','').replace(']','').lower() for header in list(blast_df['Hit Def'])]
    blast_df['Genus'] = [sciname.split(' ')[0] for sciname in list(blast_df['Scientific Name'])]
    blast_df['Species'] = [sciname.split(' ')[1] for sciname in list(blast_df['Scientific Name'])]
    blast_df['AA Sequence'] = [ aa_seq.replace('-','') for aa_seq in list(blast_df['Subject'])  ]



    ##### BIN HOMOLOG SEQUENCES. PLOT.

    bins = [50.0, 60.0, 70.0, 80.0, 90.0] 

    blast_df_50_60 = blast_df.loc[ blast_df['%Similarity'] >= 50].loc[ blast_df['%Similarity'] <60 ]
    blast_df_60_70 = blast_df.loc[ blast_df['%Similarity'] >= 60].loc[ blast_df['%Similarity'] <70 ]
    blast_df_70_80 = blast_df.loc[ blast_df['%Similarity'] >= 70].loc[ blast_df['%Similarity'] <80 ]
    blast_df_80_90 = blast_df.loc[ blast_df['%Similarity'] >= 80].loc[ blast_df['%Similarity'] <90 ]
    blast_df_90_100 = blast_df.loc[ blast_df['%Similarity'] >= 90].loc[ blast_df['%Similarity'] <100 ]

    binned_dfs = [blast_df_50_60, blast_df_60_70, blast_df_70_80, blast_df_80_90, blast_df_90_100 ]

    unique_genus_counts = [ len(b_df['Genus'].unique()) for b_df in binned_dfs]
    unique_species_counts = [ len(b_df['Species'].unique()) for b_df in binned_dfs]

    fig, [ax1, ax2, ax3] = plt.subplots(nrows=3,ncols=1, figsize= (10,5), constrained_layout = True)
    axes = fig.get_axes()

    # %Similarity hist
    ax1.hist(blast_df['%Similarity'],bins,color='mediumaquamarine')
    ax1.set_xlim(40,100)
    ax1.set_title(f'% Similarity N={len(blast_df.index)}')

    # Unique genuses bar plot
    ax2.bar(bins,unique_genus_counts,width=10,color='mediumslateblue',align='edge')
    ax2.set_xlim(40,100)
    total_genuses = len(blast_df['Genus'].unique())
    ax2.set_title(f'Unique Genuses | Total {total_genuses}')

    # Unique species bar plot
    ax3.bar(bins,unique_species_counts,width=10,color='royalblue',align='edge')
    ax3.set_xlim(40,100)
    total_species = len(blast_df['Species'].unique())
    ax3.set_title(f'Unique Species | Total {total_species}')

    fig.supxlabel("Percent Similarity Bins")
    fig.supylabel("Frequency")
    fig.savefig("Blast Summary.png",dpi=200)

    #### RANDOMLY GRAB 96x OF VARYING %Similarity. RETURN CSV
    blast_df.to_csv('blast_results.csv')

    samples = []
    for binned_df in binned_dfs: 
        
        if len(list(binned_df.index)) < 19:
            samples.append(binned_df)
            
        else:
            samples.append(binned_df.sample(19,replace=False))


    simple_cols = ['Title','%Similarity','Scientific Name','AA Sequence']
    export_df = pd.concat(samples)[ simple_cols ]
    final_df = pd.concat([export_df,pd.DataFrame({'Title':['Query'], '%Similarity': [100], 'Scientific Name': ['Query'], 'AA Sequence': [aa_sequence]})]).reset_index(drop=True)

    final_df.to_csv(f'96-well-homologs.csv')


if __name__=="__main__":
    main()