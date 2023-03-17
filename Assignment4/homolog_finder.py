import argparse
import re
import requests
import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import matplotlib as plt



'''Program to accept a protein sequence and run a BLAST query. Find homologs of the protein of different species and bin by % identity.
   Codon optimize sequences for expression in Ecoli and write to a csv.'''

ensembl_server = 'http://rest.ensembl.org'

def get_ENS_gene_id(gene,species='mouse'):
    '''Take a gene of interest and return the ensembl gene id'''       
    server = 'http://mygene.info/v3'
    query_settings = f'fields=ensembl.gene&species={species}&size=1000'
    endpoint = f'/query?q={gene}&{query_settings}'
    r = requests.get(server+endpoint)
    r = r.json()
    return r['hits'][0]['ensembl']['gene']

def get_homologs(ens_id):
    '''return dataframe of homologs of target gene'''
    #homology query
    homology_endpoint = f"/homology/id/{ens_id}?"

    r = requests.get(ensembl_server+homology_endpoint, headers={ "Content-Type" : "application/json"})
    r = r.json()
  
    #parse query
    homologous_protein_ids = [ r['data'][0]['homologies'][i]['target']['protein_id'] for i in range(0,len(r['data'][0]['homologies']))]
    perc_ids = [ r['data'][0]['homologies'][i]['target']['perc_id'] for i in range(0,len(r['data'][0]['homologies']))]
    align_seqs = [ r['data'][0]['homologies'][i]['target']['align_seq'] for i in range(0,len(r['data'][0]['homologies']))]
    perc_pos_list = [ r['data'][0]['homologies'][i]['target']['perc_pos'] for i in range(0,len(r['data'][0]['homologies']))]
    species_list = [ r['data'][0]['homologies'][i]['target']['species'] for i in range(0,len(r['data'][0]['homologies']))]

    #write df
    df = pd.DataFrame({'ProteinIDs' : homologous_protein_ids,'%Identity'  : perc_ids,
                       'Align Sequences' : align_seqs,'%Pos' : perc_pos_list, 
                       'Species' : species_list})

    return df

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
    blast_df.to_csv('blast_results.csv')

    return blast_df

def main():
    ####### ARG HANDLING #######
    parser = argparse.ArgumentParser()
    parser.add_argument("--aa_sequence", help="aa sequence of protein", type=str)
    parser.add_argument("--fasta",help="FASTA file of aa sequence of interest",type=str) #will read the last FASTA ENTRY
    parser.add_argument("--gene",help="geneID of sequence of interest",type=str)
    args = parser.parse_args()

    fasta = args.fasta
    aa_sequence = args.aa_sequence
    gene = args.gene


    if fasta:
        print(f'Reading {fasta} ...')

        # with open(fasta,'r') as fasta_file:
        #     #fasta_file.read()
        #     aa_sequence = fasta_file.read()
        # print(f'Obtained {aa_sequence} from {fasta}...')

        for seq_record in SeqIO.parse(fasta,"fasta"):
            aa_sequence = seq_record.seq
        print(f'Obtained {aa_sequence} from {fasta}...')

    elif aa_sequence:
        pass

    elif gene:
        #TODO : OBTAIN PROTEIN SEQ FROM GENE
        pass
        
    ####### BLASTP SEARCH. EXPORT AS CSV.

    blast_df = generate_blastp_df(aa_sequence)
    blast_df['%Identity'] = [ float(x/(len(aa_sequence)))*100 for x in list(blast_df['Identities'])]
    blast_df['Scientific Name'] = [re.findall('\[.*?\]',header)[-1].replace('[','').replace(']','').lower() for header in list(blast_df['Hit Def'])]
    blast_df['Genus'] = [sciname.split(' ')[0] for sciname in list(blast_df['Scientific Name'])]
    blast_df['Species'] = [sciname.split(' ')[1] for sciname in list(blast_df['Scientific Name'])]


    ##### BIN HOMOLOG SEQUENCES. PLOT.
    '''How many species within each bin? '''

    bins = [50.0, 60.0, 70.0, 80, 90]

    blast_df_50_60 = blast_df.loc[ blast_df['%Identity'] >= 50].loc[ blast_df['%Identity'] <60 ]
    blast_df_60_70 = blast_df.loc[ blast_df['%Identity'] >= 60].loc[ blast_df['%Identity'] <70 ]
    blast_df_70_80 = blast_df.loc[ blast_df['%Identity'] >= 70].loc[ blast_df['%Identity'] <80 ]
    blast_df_80_90 = blast_df.loc[ blast_df['%Identity'] >= 80].loc[ blast_df['%Identity'] <90 ]
    blast_df_90_100 = blast_df.loc[ blast_df['%Identity'] >= 90].loc[ blast_df['%Identity'] <100 ]


    binned_dfs = [blast_df_50_60, blast_df_60_70, blast_df_70_80, blast_df_80_90, blast_df_90_100 ]

    unique_genus_counts = [ len(b_df['Genus'].unique()) for b_df in binned_dfs]
    unique_species_counts = [ len(b_df['Species'].unique()) for b_df in binned_dfs]


    fig, [ax1, ax2, ax3] = plt.subplots(nrows=3,ncols=1, figsize= (10,5), constrained_layout = True)
    axes = fig.get_axes()

    # %identity hist
    ax1.hist(blast_df['%Identity'],bins,color='mediumaquamarine')
    ax1.set_xlim(42,100)
    ax1.set_title('% Identity')

    #TODO : CORRECT BAR PLOT BINNING ON XAXIS 
    # Unique genuses bar plot
    ax2.bar(bins,unique_genus_counts,width=10,color='mediumslateblue')
    ax2.set_title('Unique Genuses')


    # Unique species bar plot
    ax3.bar(bins,unique_species_counts,width=10,color='royalblue')
    ax3.set_title('Unique Species')

    fig.supxlabel("Percent Identity Bins")

    #### RANDOMLY GRAB 96x OF VARYING %identity. RETURN CSV



if __name__=="__main__":
    main()