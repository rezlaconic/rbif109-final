import requests
import json
import pandas as pd


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


####### BLASTP SEARCH





###### PARSE RESULTS




##### FILTER HOMOLOG RECORDS



##### PROCESS RECORDS AND CODON OPTIMIZE



#### RETURN CSV