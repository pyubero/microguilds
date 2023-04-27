# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:53:03 2023

@author: Pablo Yubero

This file simplifies the search of sequences for common genes (CDS and RNAs)
from NCBIs database. It takes a file with a list of organism names of interes,
and you will need to specify the query (gene name, e.g., recA) and the type,
i.e., CDS or RNA, that is because accessing NCBI is different for each.

If you have obtained manually a Refseq/Assembly pair for a specific organism,
you can add it to organismlist_correct.csv and re-run the script. It will think
that those are valid (refseq/assembly) and is it won't find the entry in the
sequences file, it will download it.

If, however, you have obtained an accession code manually, just copy/pase
the sequence to the sequences file and format nicely the header :)


It outputs three files:
    organismlist_correct.csv : It is a csv file in which each line contains the
    sufficient information to recover the source of the sequence obtained.

    organismlist_failed.csv : It is a list of organisms that have failed so you
    can easily identify mispellings for example.

    sequences_QUERY.fasta : A fasta file containing the sequences obtained, and
    a simplified header with source information.

For example you can run the script as:

python get_sequences_from_ncbi.py organislist.csv -q recA -t CDS --verbose --clear_all
"""
import os
import wget
import requests
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from bs4 import BeautifulSoup
from time import sleep


def name_to_refseq(organism_name : str, gene = "genome"):
    BASE_URL = 'https://www.ncbi.nlm.nih.gov/nuccore/?term='

    # Create search term
    search_term = organism_name.replace(' ','+')+'+'+gene

    # Making a GET request
    sleep(1/3)
    r = requests.get( BASE_URL+search_term )

    # Parse request
    soup = BeautifulSoup(r.content, 'html.parser')

    # Find object that contains relevant info
    s = soup.find('ul', class_='ncbi-inline-list')
    if (s is None) or ('RefSeq' not in s.text):
        return None

    # Get your RefSeq
    return s.li.text.split(' ')[1]


def refseq_to_assembly( refseq ):
    BASE_URL = 'https://www.ncbi.nlm.nih.gov/assembly/REFSEQ/'

    if refseq is None:
        return None, None

    # Making a GET request
    sleep(1/3)
    r = requests.get( BASE_URL.replace('REFSEQ', refseq))

    # Parse request
    soup = BeautifulSoup(r.content, 'html.parser')

    # Find object that contains relevant info
    s = soup.find('h1', class_="marginb0 margin_t0")
    if (s is None):
        return None, None
    assembly = s.getText()

    # Extract real name of assembly's organism
    s = soup.find(id="summary").find_all('dd')
    for _ in s:
        if "Organism name" in _.previous_sibling.get_text(strip=True):
            name = ' '.join( _.get_text(strip=True).split(' ')[:2] )

    return assembly, name


def download_assembly_from_ncbi(refseq, assembly, type_, filename_out=None, bar=True):
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
    filename = refseq+'_'+assembly_name+'_%s_from_genomic.fna.gz' % type_.lower()
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'+post+'_'+assembly_name+'/'+filename

    if bar:
        bar = bar_custom
    else:
        bar = bar_none

    if os.path.exists(filename_out):
        os.remove(filename_out) # if exist, remove it directly

    try:
        response = wget.download(url, filename_out, bar=bar_custom)
        return filename_out
    except:
        return False



def sequence_from_assembly(query, refseq, assembly, type_):

    if (refseq is None) or (assembly is None):
        return '',''

    filename_temp = "REFSEQ_data.temp".replace('REFSEQ', refseq )
    filename_assembly = './TYPE/REFSEQ.fna.gz'.replace("TYPE", type_).replace('REFSEQ', refseq)

    create_dir_if_doesnt_exist("./TYPE".replace("TYPE", type_) )


    # Download assembly genomes, and look for the sequences of interest
    if not download_assembly_from_ncbi(refseq, assembly, type_, filename_out=filename_assembly):
        if VERBOSE:
            print('<W> Could not download the assembly code of %s' % ( refseq ) )
        return '' , ''


    # Unzip file
    gunzip(filename_assembly, filename_temp)

    # Parse file searching for the query name
    header, body = get_sequence_from_fasta( filename_temp,  "%s" % query)

    if os.path.exists(filename_temp):
        os.remove(filename_temp)

    return header, body


def get_sequence_from_fasta( filename : str, query : str):
    header, body = '', ''
    with open(filename,'r') as file:
        line = 'X'
        while line != '':
            line = file.readline()

            if query in line:
                header = line.replace('\n','')
                body =''
                line = file.readline()
                while (len(line)>0) and (line[0]!='>'):
                    body += line
                    line = file.readline()
                break
    return header, body


def format_header( header, query,  name, truename, refseq, assembly, accession ):
    loc = header_get_location(header)
    if truename is None:
        truename = name
    header = ">{}:{} [original_query={}] [refseq={}] [assembly={}] [seq_query={}] [loc={}] [accession={}]\n".format(
        truename.replace(' ','_'), query, name, refseq, assembly, query, loc, accession)
    return  header


def header_get_location( header : str):
    if "location" in header:
        return header.split('location')[1].split("]")[0][1:]
    else:
        return "Unknown"


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


def create_dir_if_doesnt_exist( path ):
    if os.path.isfile(path):
        raise PathError("Not a dir path.")

    if not os.path.isdir(path):
        os.mkdir(path)


def sequence_from_entrez_search( query ):
    idlist = NCBI_search( query )
    fastas = NCBI_fetch( idlist[:5] )

    if fastas is None:
        return '', ''

    for fasta in fastas:
        _h, _b = fasta
        if every_word_in_str( query.split(' '), _h):
            header = _h
            body = _b
            return header, body
    return '', ''

def every_word_in_str( query : list, reference : str):
    return np.all( [word in reference for word in query] )


def NCBI_search(term, db='nuccore', retmode='json'):
    URL_ESEARCH='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    PARAMS_SEARCH = {
        'db' : db, # or nucleotide
        'term' : term,
        'retmode' : retmode,
        'sort' : 'relevance' }

    sleep(1/3)
    response = requests.get( URL_ESEARCH, params = PARAMS_SEARCH, timeout = 10)

    if response.status_code==200:
        json = response.json()
        idlist = json['esearchresult']['idlist']
        return idlist
    else:
        print('<W> Error in the NCBI search.')
        return None


def NCBI_fetch(GI, db='nuccore', rettype='fasta',retmode='text'):
    if GI is None:
        return None

    URL_EFETCH= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    PARAMS_FETCH = {
        'db' : db, # or nucleotide
        'id' : GI,
        'rettype' : rettype,
        'retmode' : retmode}

    sleep(1/3)
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


def add_data(data, name=None, truename=None, refseq=None, assembly=None, accession=None, success=[]):
    data["Name"].append( name )
    data["TrueName"].append( truename )
    data["Refseq"].append( refseq )
    data["Assembly"].append( assembly )
    data["Accession"].append( accession )
    data["SuccessWith"].append( success )
    return data



def is_sequence_in_fasta(filename, query_list):
    with open(filename,'r') as file:
        for line in file.readlines():
            if every_word_in_str( query_list , line):
                return True
    return False


def export_failed(organism, filename):
    with open(filename,'a+') as file:
        file.write(organism+'\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        prog = 'get_sequences_from_ncbi.py',
                        description = "get_sequences_from_ncbi.py belongs to the microguilds package. It searches for Refseqs and assembly codes in NCBI's database from organism names.",
                        epilog = 'Please visit github.com/pyubero/microguilds for more details.')

    parser.add_argument('input', metavar='filename')
    parser.add_argument('-q','--query', metavar='sequence_name', help="Name of genes or sequences queried.")
    parser.add_argument('-t','--type', metavar="type_", help="Either CDS or RNA", choices=["CDS", "RNA"])
    parser.add_argument('-v', '--verbose', action='store_true', help="Inline print out of detailed progress.")  # on/off flag
    parser.add_argument('--clear_all', action='store_true', help='Clears all output files before running the script.')

    args = parser.parse_args("organislist.csv -q potF -t CDS --verbose --clear_all")


    # External variables
    FILENAME_IN = args.input
    QUERY = args.query
    TYPE = args.type
    CLEAR_ALL = args.clear_all
    VERBOSE = args.verbose

    # Internal variables
    _nerrors = 0
    _out_table ='organismlist_correct.csv'
    _out_failed='organismlist_failed.csv'
    _out_fasta ='sequences_{}.fasta'.format(QUERY)
    _data = {"Name" : [],
             "TrueName" : [],
             "Refseq" : [],
             "Assembly" : [],
             "Accession" : [],
             "SuccessWith" : []}

    # Decrlate verbosity
    if VERBOSE:
        def verboseprint(*args, **kwargs):
            print(*args, **kwargs)
    else:
        verboseprint = lambda *a, **k: None # do-nothing function


    if CLEAR_ALL:
        if os.path.exists(_out_table):
            os.remove(_out_table)
        if os.path.exists(_out_fasta):
            os.remove(_out_fasta)
    if os.path.exists(_out_failed):
        os.remove(_out_failed)


    # 1. Load organism names
    orgnames = pd.read_csv(FILENAME_IN , header=0)["Name"]
    verboseprint("Found %d input organisms." % len(orgnames))

    #... load previous session
    if os.path.exists(_out_table):
        _data = pd.read_csv(_out_table,keep_default_na=False).to_dict()
        _data["Name"] = [value for key, value in _data["Name"].items() ]
        _data["TrueName"] = [value for key, value in _data["TrueName"].items() ]
        _data["Refseq"] = [value for key, value in _data["Refseq"].items() ]
        _data["Assembly"] = [value for key, value in _data["Assembly"].items() ]
        _data["Accession"] = [value for key, value in _data["Accession"].items() ]
        _data["SuccessWith"] = [value for key, value in _data["SuccessWith"].items() ]


    # ... for every organism in the list
    for jj, organism in enumerate(orgnames):

        verboseprint("--",jj,"-- ", organism)
        header, body = "", ""
        truename, refseq, assembly, accession = "","","",""

        # Resume if already computed
        if organism in _data["Name"]:
            idx = np.argwhere(organism == np.array(_data["Name"])).flatten()[0]
            truename, refseq, assembly, accession= [ _data[key][idx] for key in _data ][1:-1]

            if is_sequence_in_fasta(_out_fasta, [organism, refseq, assembly, accession]):
                verboseprint("Already found in ", _out_table)
                verboseprint("")
                continue

        # 2. Obtain Refseq from name
        if refseq == "":
            refseq = name_to_refseq(organism)
        verboseprint('Refseq: ', refseq)


        # 3. Obtain assembly from refseq
        if assembly == "":
            assembly, truename = refseq_to_assembly( refseq )
        verboseprint('Assembly: ', assembly)
        verboseprint('Tru name: ', truename)


        # 4a. Search sequence(s) in assembly genome
        if (refseq !="") and (assembly !=""):
            header, body = sequence_from_assembly( QUERY, refseq , assembly, TYPE )
        verboseprint("Seq length: ", len(body.replace('\n','')))


        # 4b. If empty, search in NCBI's nuccore
        if body=="":
            header, body = sequence_from_entrez_search( organism+' '+QUERY )
            accession = header.split(' ')[0][1:]
            verboseprint("Using EUTILS...")
            verboseprint("Accession: ", accession)
            verboseprint("Seq length: ", len(body.replace('\n','')))


        # 5a. If failed...
        if body=="":
            _nerrors += 1
            verboseprint('')
            export_failed(organism, _out_failed)
            continue


        # 5b. If successful...
        header = format_header( header, QUERY,  organism, truename, refseq, assembly, accession )
        with open(_out_fasta,'a+') as file:
            file.write(header)
            file.write(body)


        # 6. Update output table
        _data = add_data(
            _data,
            name = organism,
            truename= truename,
            refseq = refseq,
            assembly=assembly,
            accession=accession,
            success=[]        )

        pd.DataFrame.from_dict(_data).to_csv(_out_table, index=False)

        verboseprint('')

verboseprint("Found {} errors.".format(_nerrors) )
