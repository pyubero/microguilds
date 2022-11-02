# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 17:48:09 2022

@author: Pablo

# This code proceeds as follows:
    1. Using an organism name, it searches for the RefSeq code of its assembly genome
    2. It searches for the latest version of the assembly genome
    3. It downloads the assembly genome in CDS format
    [nope] 4. It unzips the file
    [nope] 5. It parses the file, searches for the gene of interest and saves it
"""
import os
import wget
import requests
from tqdm import tqdm
from bs4 import BeautifulSoup

def print_to_log( data, flog='log.log' ):
    with open(flog , 'a+') as file:
        file.write(data+'\n')
    print(data)

def name2refSeq(organism_name):
    BASE_URL = 'https://www.ncbi.nlm.nih.gov/nuccore/?term='


    # Create search term
    search_term = orgname.replace(' ','+')+'+genome'
    
    # Making a GET request
    r = requests.get( BASE_URL+search_term )
    
    # Parse request
    soup = BeautifulSoup(r.content, 'html.parser')
    
    # Find object that contains relevant info
    s = soup.find('ul', class_='ncbi-inline-list')
    if (s is None) or ('RefSeq' not in s.text):
        # print('Could not find anything for %s' % orgname )
        return None
    
    # Get your RefSeq
    answer = s.li.text.split(' ')[1]    
    return answer


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




def is_str_in_file(filename, string):
    with open(filename,'r') as file:
        for line in file.readlines():
            if string in line:
                return True
    return False


FILENAME_IN ='organismpotf.csv'
# FILENAME_OUT='organisms_data.txt'
FILENAME_OUT='data_recA.csv'
FILENAME_LOG = 'refseq_scraping.log'

RESUME = True
#...
_nerrors = 0




# Get list of organisms' names
organisms_names = []
with open(FILENAME_IN) as file:
    for line in file.readlines():
        organisms_names.append( line.replace('\n','') )



# Clear output files
if not RESUME:
    with open(FILENAME_OUT , 'w+') as file:    pass
    with open('log.log' , 'w+') as file:    pass


# Search and download RefSeqs
for ii, orgname in tqdm( enumerate(organisms_names) , total=len(organisms_names)) :

    if RESUME and is_str_in_file(FILENAME_OUT, orgname):
        continue

    # Step 1: Find RefSeq code of organism from the organism name
    refseq = name2refSeq(orgname)
    # ... if it can not find any, save anyways with empty cells
    if refseq is None:
        _nerrors += 1
        print_to_log('<W> Could not find the RefSeq of %s' % orgname,
                     flog = FILENAME_LOG)
    
        # Create organism data file
        with open(FILENAME_OUT,'a+') as file:
            file.write('%s,,\n' % (orgname,) ) 

        continue
    
   
    
    # Step 2: Find assembly version and code from the refseq
    assembly_name = refseq2assembly( refseq )
    # ... if it can not find any, save anyways with empty cell
    if assembly_name is None:
        _nerrors += 1
        print_to_log('<W> Could not find the assembly code of %s, %s' % (orgname,refseq) ,
                     flog = FILENAME_LOG)
        # Create organism data file
        with open(FILENAME_OUT,'a+') as file:
            file.write('%s, %s,\n' % (orgname, refseq) ) 
                
        continue


    # Append to organism data file if everything went OK
    with open(FILENAME_OUT,'a+') as file:
        file.write('%s, %s, %s\n' % (orgname, refseq, assembly_name) ) 
            
     
    
    


print_to_log('Finished!', flog = FILENAME_LOG)
print_to_log('Had %d errors' % _nerrors, flog = FILENAME_LOG)





