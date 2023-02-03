# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:11:04 2022

@author: Pablo


# Try importing organisms_data using a pandas dataframe
"""

import os
import wget
import requests
import argparse
from tqdm import tqdm
from bs4 import BeautifulSoup
from time import sleep




def load_refseqs_and_assemblies_from_file(filename : str):

	# Get list of organisms' names
	names = []
	refseqs = []
	assemblies = []
	with open(filename) as file:
	    for line in file.readlines():
	        line = line.replace('\n','')
	        _name   = line.split(',')[0]
	        _refseq = line.split(',')[1].replace(' ','')
	        _assemb = line.split(',')[2].replace(' ','')

	        if len(_assemb) > 1:
	            names.append( _name )
	            refseqs.append( _refseq)
	            assemblies.append( _assemb)
	        else:
	        	if VERBOSE: print("Unable to fetch assembly code of %s" % _name )

	return names, refseqs, assemblies      




def download_assembly_from_ncbi(refseq, assembly, type, filename_out=None, bar=True):
    def bar_custom(current, total, width=80):
        prefix = '\b'*999
        print(prefix + "Downloadong: %d%% [%d / %d] bytes" % (current / total * 100, current, total), end='')
    def bar_none(current, total):
        pass


	if filename_out is None:
        filename_out = "%s.fna.gz" % refseq

    if os.path.exists(filename_out):
    	if VERBOSE: print('File %s already exists.' % filename_out)
    	return True

    if assembly is None:
        assembly_name = refseq2assembly(refseq)
    else:
        assembly_name = assembly
        
    assembly_name = assembly_name.replace(' ','_')
    header = refseq.split('_')[0]
    s1 = refseq.split('_')[1][0:3]
    s2 = refseq.split('_')[1][3:6]
    s3 = refseq.split('_')[1][6:9]
    
    post = '/'.join([header,s1,s2,s3,refseq]) 
    filename = refseq+'_'+assembly_name+'_%s_from_genomic.fna.gz' % type.lower()    
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'+post+'_'+assembly_name+'/'+filename
    
    if bar:  bar = bar_custom
    else:  bar = bar_none
    
    if os.path.exists(filename_out):
        os.remove(filename_out) # if exist, remove it directly
        
    try:
        response = wget.download(url, filename_out, bar=bar_custom) 
        return True
    except:
        return False


def sequence_from_assembly( filename : str, query : str):
    header, body = '', ''
    with open(filename,'r') as file:
        line = 'X'
        while line!='':
            line = file.readline()
            
            if query in line:
                header = line.replace('\n','') + " [orgname=%s] [refseq=%s]\n" % (_name,_refseq)
                body =''
                line = file.readline()
                while (len(line)>0) and (line[0]!='>'):
                    body += line
                    line = file.readline()
    
                break
    return header, body

     
def every_word_in_str( query : list, reference : str):
	return np.all( [word in reference for word in query] )


def sequence_from_entrez_search( query ):
	idlist = NCBI_search( query )
    fastas = NCBI_fetch( idlist[:5] )            
    for fasta in fastas:
        _h, _b = fasta
        if every_word_in_str( query.split(' '), _h):
            header = _h+'\n'
            body = _b
            return header, body


def NCBI_search(term, db='nuccore', retmode='json'):
    URL_ESEARCH='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    PARAMS_SEARCH = {
        'db' : db, # or nucleotide
        'term' : term,
        'retmode' : retmode,
        'sort' : 'relevance' }
    
    response = requests.get( URL_ESEARCH, params = PARAMS_SEARCH, timeout = 10)
    if response.status_code==200:
        json = response.json()
        idlist = json['esearchresult']['idlist']
        sleep(0.2)
        return idlist
    else:
        print('<W> Error in the NCBI search.')
        return None


def NCBI_fetch(GI, db='nuccore', rettype='fasta',retmode='text'):
    URL_EFETCH= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    PARAMS_FETCH = {
        'db' : db, # or nucleotide
        'id' : GI,
        'rettype' : rettype,
        'retmode' : retmode}
    response = requests.get( URL_EFETCH, params = PARAMS_FETCH, timeout = 5)
    if response.status_code==200:
        text = response.text.split('>')[1:]
        fastas = []
        for mixed in text:
            lines = mixed.split('\n')
            header = '>'+lines[0]
            body = '\n'.join(lines[1:])
            fastas.append( [header, body] )
        
        return fastas


#if __name__ == '__main__':

parser = argparse.ArgumentParser(
                    prog = 'fetch_sequences.py',
                    description = 'This program retrieves the nucleotide sequences of a given gene or RNA given a list of organism names, refseqs and assembly codes. Typically used in combination with other microbial_guilds scripts.',
                    epilog = 'Please visit the github repo for further details.')


parser.add_argument('filename', help="Input file with a list of organism names, refseqs and assembly codes separated by commas.")           # positional argument
parser.add_argument('qname', help="Query name of the cds or rna sequence.")      # option that takes a value
parser.add_argument('-type', help="Sequence type, either CDS or RNA", choices=["CDS", "RNA"])
parser.add_argument('-min', help="Minimum length of sequences. It can really improve results to set a reasonable limit.")
parser.add_argument('--clear_all', action='store_true', help="Clears output files before running the script")
parser.add_argument('--keep_temp', action='store_true', help="Prevents deleting temporal CDS files.")
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

args = parser.parse_args("organisms_data.csv potF -type CDS -min 1200 --keep_temp --verbose".split(' '))


# Input variables #
QNAME = args.qname
TYPE  = args.type
FILENAME_DATA = args.filename
MINLEN = args.min 
CLEAR_ALL = args.clear_all
KEEP_TEMP = args.keep_temp 
VERBOSE = args.verbose 

# Internal variables #

FILENAME_OUT = 'seq_%s.fasta' % QNAME
FILENAME_TEMP = "REFSEQ_data.temp"
FILENAME_ASSEMBLY = './TYPE/REFSEQ.fna.gz'.replace("TYPE", TYPE)
#...
_nerrors = 0
        


if CLEAR_ALL:
    with open(FILENAME_OUT , 'w+') as file:    pass


names, rqfseqs, assemblies = load_refseqs_and_assemblies_from_file(FILENAME_DATA)
if VERBOSE: print("Loaded names, refseqs and assembly codes for %d organisms." % len(names) )


for ii in range(len(names)):
	
	# Download assembly genomes, and look for the sequences of interest        
	if not download_assembly_from_ncbi(refseqs[ii], assemblies[ii], TYPE, filename_out=None):
	    _nerrors += 1
        if VERBOSE: print('<W> Could not download the assembly code of %s, %s' % (names[ii], refseqs[ii] ) )
	    continue
    

    # Unzip file
    filename_assembly = FILENAME_ASSEMBLY.replace('REFSEQ', refseqs[ii])
    filename_temp = FILENAME_TEMP.replace('REFSEQ', refseqs[ii] )
    gunzip(filename_assembly, filename_temp)
    
    # Parse the genome file and keep what is important
    # Parse file searching for the query name
    header, body = sequence_from_assembly( filename_temp,  "%s" % QNAME)

    if not KEEP_TEMP:
        if os.path.exists(filename_assembly):
            os.remove(filename_assembly) 
        if os.path.exists(filename_temp):
            os.remove(filename_temp) 


    # If the traditional way fails, lets try some more...            
    if (len(header) < 2) or (len(body) < MINLEN):
    	header, body = sequence_from_entrez_search( names[ii]+' '+QNAME )

    # If it was succesfull
    if len(header)>1 and len(body)>MINLEN:
        if VERBOSE: print('[%3d/%3d] %s, %s, %s, %d' % (ii, len(names), names[ii], refseqs[ii], assemblies[ii], len(body.replace('\n','') )) )
        
        with open(FILENAME_OUT,'a+') as file:
            file.write(header)
            file.write(body)  
        
    else:
        _nerrors+=1
        print('<W> Could not find the %s of %s, %s' % (QNAME, names[ii], refseqs[ii]))
            

print_to_log('Finished!')
print_to_log('Had %d errors' % _nerrors)