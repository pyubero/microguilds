# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:11:04 2022

@author: Pablo
"""

import os
import wget
import requests
from tqdm import tqdm
from bs4 import BeautifulSoup
from time import sleep

def print_to_log( data, flog='log.log' ):
    with open(flog , 'a+') as file:
        file.write(data+'\n')
    print(data)


def downloadAssembly(refseq, assembly=None, filename_out=None, bar=True):
    def bar_custom(current, total, width=80):
        prefix = '\b'*999
        print(prefix + "Downloadong: %d%% [%d / %d] bytes" % (current / total * 100, current, total), end='')
    def bar_none(current, total):
        pass

    if assembly is None:
        assembly_name = refseq2assembly(refseq)
    else:
        assembly_name = assembly
        
    if filename_out is None:
        filename_out = "%s.fna.gz" % refseq
        
    assembly_name = assembly_name.replace(' ','_')
    header = refseq.split('_')[0]
    s1 = refseq.split('_')[1][0:3]
    s2 = refseq.split('_')[1][3:6]
    s3 = refseq.split('_')[1][6:9]
    
    post = '/'.join([header,s1,s2,s3,refseq]) 
    # filename = refseq+'_'+assembly_name+'_cds_from_genomic.fna.gz'
    filename = refseq+'_'+assembly_name+'_rna_from_genomic.fna.gz'
    
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'+post+'_'+assembly_name+'/'+filename
    
    if bar:
        bar = bar_custom
    else:
        bar = bar_none
    
    if os.path.exists(filename_out):
        os.remove(filename_out) # if exist, remove it directly
        
    try:
        response = wget.download(url, filename_out, bar=bar_custom) 
        return True
    except:
        return None
    
def refseq2assembly( refseq):
    BASE_URL = 'https://www.ncbi.nlm.nih.gov/assembly/REFSEQ/'
    
    # Making a GET request
    r = requests.get( BASE_URL.replace('REFSEQ', refseq))
    
    # Parse request
    soup = BeautifulSoup(r.content, 'html.parser')
    
    # Find object that contains relevant info
    s = soup.find('h1', class_="marginb0 margin_t0")
    if (s is None):
        return None
    
    # Get your RefSeq
    return s.getText()    




def gunzip(source_filepath, dest_filepath, block_size=65536):
    import gzip

    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)
        
def is_str_in_file(filename, string):
    if os.path.exists(filename):
        with open(filename,'r') as file:
            for line in file.readlines():
                if string in line:
                    return True
    return False
    

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
        
        


RNA = '16S'
MINLEN = 1200
# FILENAME_DATA='organisms_data.txt'
FILENAME_DATA='data_recA.csv'
FILENAME_OUT = 'seq_%s.fasta' % RNA
FILENAME_LOG = 'download_rna.log'
REMOVE_GZIP = False
REMOVE_CDS  = True
RESUME = True
#...
FILENAME_TEMP = "REFSEQ_rna.temp"
FILENAME_ASSEMBLY = './RNA/REFSEQ.fna.gz'
#...
_nerrors = 0
        


if not RESUME:
    with open(FILENAME_OUT , 'w+') as file:    pass
    with open(FILENAME_LOG, 'w+') as file:    pass




# Get list of organisms' names
names = []
refseqs = []
assemblies = []
with open(FILENAME_DATA) as file:
    for line in file.readlines():
        line = line.replace('\n','')
        _name   = line.split(',')[0]
        _refseq = line.split(',')[1].replace(' ','')
        _assemb = line.split(',')[2].replace(' ','')

        if len(_assemb) > 1:
            names.append( _name )
            refseqs.append( _refseq)
            assemblies.append( _assemb)
        
        
        
# Download assembly genomes, and look for the sequences of interest        
for ii in range(len(names)):
    _name = names[ii]
    _refseq= refseqs[ii]
    _assemb= assemblies[ii]
    
    if RESUME and is_str_in_file(FILENAME_OUT, _name):
        continue
    
        
    # Step 3: Download assembly genome and the CDS file in fasta
    filename_assembly = FILENAME_ASSEMBLY.replace('REFSEQ', _refseq)
    
    suc = True
    if not os.path.exists(filename_assembly):
        suc= downloadAssembly( _refseq, 
                               assembly = _assemb, 
                               filename_out = filename_assembly,
                               bar = False)
        sleep(1)
    # print('')
    
    if not suc:
        _nerrors += 1
        print_to_log('<W> Could not download the assembly code of %s, %s' % (_name, _refseq),
                     flog=FILENAME_LOG )
        continue
    
    # Step 4: Unzip file
    filename_temp = FILENAME_TEMP.replace('REFSEQ', _refseq)
    gunzip(filename_assembly, filename_temp)
    
    # Step 5: Parse the genome file and keep what is important
    # Parse file searching for the RNA
    search_term = "%s" % RNA
    header, body = '', ''
    
    with open(filename_temp,'r') as file:
        line = 'X'
        while line!='':
            line = file.readline()
            
            if search_term in line:
                header = line.replace('\n','') + " [orgname=%s] [refseq=%s]\n" % (_name,_refseq)
                body =''
                line = file.readline()
                while (len(line)>0) and (line[0]!='>'):
                    body += line
                    line = file.readline()
    
                break

    
    if REMOVE_GZIP:
        if os.path.exists(filename_assembly):
            os.remove(filename_assembly) 
        
    if REMOVE_CDS:
        if os.path.exists(filename_temp):
            os.remove(filename_temp) 
     
            
    # If the traditional way fails, lets try some more...            
    if len(header)<2 or len(body)<MINLEN:
        idlist = NCBI_search( _name+' '+RNA )
        fastas = NCBI_fetch( idlist[:5] )            
        for fasta in fastas:
            _h, _b = fasta
            if (_name.split(' ')[0] in _h) and (_name.split(' ')[1] in _h) and (RNA in _h):
                header = _h+'\n'
                body = _b
                break
            
            

            
    if len(header)>1 and len(body)>MINLEN:
        print_to_log('[%3d/%3d] %s, %s, %s, %d' % (ii, len(names), _name, _refseq, _assemb, len(body.replace('\n','') )),
                      flog=FILENAME_LOG )
        
        with open(FILENAME_OUT,'a+') as file:
            file.write(header)
            file.write(body)  
        pass
    else:
        _nerrors+=1
        print('<W> Could not find the 16S of %s, %s' % (_name, _refseq))
            
    # print('')



print_to_log('Finished!')
print_to_log('Had %d errors' % _nerrors)




