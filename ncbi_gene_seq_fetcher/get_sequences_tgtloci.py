# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 10:06:47 2022

@author: Pablo Yubero

This script is to be used with files downloaded from the TargetedLoci db of ncbi/ftp:
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/

For example, call it as:
    python get_sequences_tgtloci.py organisms_input.csv bacteria.16SrRNA.fna 16S -v
"""

import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

def get_sequence(search_terms, fastafile):
    
    # First check that search term is provided as a list
    assert type([search_terms,])==list
    
    with open(fastafile,'r') as file:
        headers = []
        bodies  = []
        header = ''
        body = ''
        _saving = False
        
        for line in file.readlines():
            if line[0]!='>' and _saving:
                body += line
            
            if line[0]=='>':
                # Here we have a fully working (header, body) pair
                if np.all([_ in header for _ in search_terms]):
                    headers.append(header)
                    bodies.append(body)
                
                #... renew header and body
                header = line
                body = ''
                _saving = True
                
        return headers, bodies
        

def format_header( header, query,  name, truename, refseq, assembly, accession ):
    if truename is None:
        truename = name
    header = ">{}:{} [original_query={}] [refseq={}] [assembly={}] [seq_query={}] [loc={}] [accession={}]\n".format(
        truename.replace(' ','_'), query, name, refseq, assembly, query, "", accession)
    return  header


def define_translations():
    name_translations = {}
    name_translations.update( {'Ruegeria litorea'           : 'Falsiruegeria litorea'} )
    name_translations.update( {'Ruegeria pomeyori'          : 'Ruegeria pomeroyi'} )
    name_translations.update( {'Tritonibacter mobile'       : 'Tritonibacter mobilis'} )
    name_translations.update( {'Aliiseovarius crassoterae'  : 'Aliiroseovarius crassostreae'} )
    name_translations.update( {'Zobellela maritima'         : 'Zobellella maritima'} )
    name_translations.update( {'Salipiger pacificus'        : 'Yangia pacifica'} )
    name_translations.update( {'Halocynthiibacter arcticus' : 'Falsihalocynthiibacter arcticus'} )
    name_translations.update( {'Labrenzia alba'             : 'Roseibium album'} )
    name_translations.update( {'Mesorhizobium oceanicum'    : 'Aquibium oceanicum'} )
    name_translations.update( {'Pseudomonas oceani'         : 'Halopseudomonas oceani '} )
    name_translations.update( {'Pseudomonas aestusnigri'    : 'Halopseudomonas aestusnigri '} )
    name_translations.update( {'Pseudopelagicola gijangensis'   : 'Shimia gijangensis'} )
    name_translations.update( {'Thalassobius aestuarii'         : 'Shimia aestuarii'} )
    name_translations.update( {'Roseivivax pacificus'           : 'Allosediminivita pacifica'} )
    name_translations.update( {'Puniceibacterium sediminis'     : 'Pseudopuniceibacterium sediminis'} )
    name_translations.update( {'Litorimicrobium taeanense'      : 'Thalassobius taeanensis'} )
    name_translations.update( {'Pseudoruegeria haliotis'        : 'Aliiruegeria haliotis'} )
    name_translations.update( {'Tropicibacter phthalicicus'     : 'Pelagimonas phthalicica'} )
    name_translations.update( {'Ruegeria mobilis'               : 'Tritonibacter mobilis'} )
    name_translations.update( {'Labrenzia aggregata'            : 'Roseibium aggregatum'} )
    name_translations.update( {'Labrenzia alexandrii'           : 'Roseibium alexandrii'} )
    name_translations.update( {'Nesiotobacter exalbescens'      : 'Pseudovibrio exalbescens'} )
    name_translations.update( {'Methylobacterium salsuginis'    : 'Methylorubrum salsuginis'} )
    name_translations.update( {'Desulfobulbus japonicus'        : 'Desulfogranum japonicum'} )
    name_translations.update( {'Pseudomonas pachastrellae'      : 'Halopseudomonas pachastrellae'} )
    name_translations.update( {'Pseudomonas pelagia'            : 'Halopseudomonas pelagia'} )
    name_translations.update( {'Balneatrix sp003194385'         : 'Balneatrix'} )
    name_translations.update( {'Pseudomonas zhaodongensis'      : '[Pseudomonas] zhaodongensis'} )
    name_translations.update( {'Amylibacter sp900197625'        : 'Amylibacter'} )
    name_translations.update( {'Anderseniella sp900149695'      : 'Anderseniella'} )
    name_translations.update( {'Pseudoruegeria lutimaris'       : 'Aliiruegeria lutimaris'} )
    name_translations.update( {'Azoarcus toluclasticus'         : 'Aromatoleum toluclasticum'} )
    name_translations.update( {'Janthinobacterium sp002127585'  : 'Janthinobacterium'} )
    # name_translations.update( {'Cognatiyoonia sediminum'        : ''} )
    # name_translations.update( {'Pseudorhizobium pelagicum'      : ''} )
    # name_translations.update( {'Motiliproteus coralliicola'     : ''} )
    # name_translations.update( {'Bermanella sp002683575'         : 'Bermanella'} )
    # name_translations.update( {'Zobellella maritima'            : ''} )
    # name_translations.update( {'Enterobacter himalayensis'      : ''} )# name_translations.update( {'Grimontia indica'               : ''} )
    return name_translations


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog = 'get_sequences_tgtloci.py',
                        description = "get_sequences_16S.py belongs to the microguilds package. It searches for sequences from a TargetedLoci file (NCBIs database) of a list of organisms.",
                        epilog = 'Please visit github.com/pyubero/microguilds for more details.')
    
    parser.add_argument('input', metavar='filename', help="Input file with the list of organisms of interest.") 
    parser.add_argument('sequencesfile', metavar='sequencesfile', help="Path to file with all reference sequences.") 
    parser.add_argument('name', metavar="name", help="Name of sequence queried, necessary for the formatting of the output fasta.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Inline print out of detailed progress.")  # on/off flag

    
    args = parser.parse_args( )
    # args = parser.parse_args( "organisms_input.csv bacteria.16SrRNA.fna 16S -v".split(' '))
    
    
    # External variables #
    FILENAME_IN = args.input
    FASTAFILE = args.sequencesfile
    QNAME = args.name
    FILENAME_OUT='sequences_16s.fasta'
    VERBOSE = True
    CLEAR_ALL=True

    # Internal variables # 
    _translations = define_translations()


    # Decrlate verbosity
    if VERBOSE:
        def verboseprint(*args, **kwargs):
            print(*args, **kwargs)
    else:
        verboseprint = lambda *a, **k: None # do-nothing function


    # Clear all    
    if CLEAR_ALL:
        with open(FILENAME_OUT , 'w+') as file:    pass    


    # 1. Load organism names
    orgnames = pd.read_csv(FILENAME_IN , header=0)["Name"]
    verboseprint("Found %d input organisms." % len(orgnames))
    
    
    for name in tqdm(orgnames):
        
        # 2. Obtain query name
        if name in _translations:
            search_term = _translations[name]
        else:
            search_term = name
    
        # 3. Get sequence
        headers, bodies = get_sequence( [search_term,] , FASTAFILE )
    
        if len(headers)==0:
            verboseprint('\n<W> Did not find a {} for {}'.format(QNAME, name) )
            continue
        elif len(headers)>=1:
            lengths = np.array( [ len(_) for _ in bodies ])
            idx = np.argmax(lengths)
        else:
            idx = 0
    
        # Retrieve header
        header, body = headers[idx], bodies[idx]
        
        # Format header
        truename = ' '.join( header.split(' ')[1:3])
        accession = header.split(' ')[0][1:]
        header = format_header( None, QNAME,  name, truename, "", "", accession )
        
        # Export headers[idx] and bodies[idx]
        with open(FILENAME_OUT,'a+') as file:
            #...
            file.write( header)
            file.write( body )
    
    
    
    
    








