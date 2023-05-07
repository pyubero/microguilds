# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:53:03 2023

@author: Pablo Yubero

This file simplifies the search of sequences for common genes (CDS and RNAs)
from NCBIs database. It takes a file with a list of organism names of interest,
and a query (gene name, e.g., recA) and its type (CDS or RNA), that is because
access to NCBI is different in each case.

It outputs three files:
    - organismlist_correct.csv : It is a csv file in which each line contains
    the sufficient information to recover the source of the sequence obtained.

    - organismlist_failed.csv : It is a list of organisms that have failed so
    you can easily identify mispellings for example.

    - QUERY_sequences.fasta : A fasta file containing the sequences obtained,
    and a simplified header with source information.

For some failed organisms, you may manually obtain a Refseq/Assembly pair,
you can add them to organismlist_correct.csv and re-run the script. It will
identify those as valid (refseq/assembly) and as it won't find the entry in the
sequences file, it will automatically retrieve them.

If, however, you have obtained a sequence manually, just copy/paste the
sequence to the sequences file and format the header nicely :)

** Example **
python get_sequences_from_ncbi.py organislist.csv -q recA -t CDS --verbose
python get_sequences_from_ncbi.py organislist.csv -q 16S -t RNA --verbose
"""
import os
import argparse
from time import sleep
import wget
import requests
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup


def name_to_refseq(organism_name: str, gene="genome"):
    '''Queries NCBI's nuccore to fetch the refseq given a species name.'''
    base_url = 'https://www.ncbi.nlm.nih.gov/nuccore/?term='

    # Create search term
    search_term = organism_name.replace(' ', '+') + '+' + gene

    # Making a GET request
    sleep(1/2)
    get_request = requests.get(base_url + search_term, timeout=10)

    # Parse request
    soup = BeautifulSoup(get_request.content, 'html.parser')

    # Find object that contains relevant info
    soup_ul = soup.find('ul', class_='ncbi-inline-list')
    if (soup_ul is None) or ('RefSeq' not in soup_ul.text):
        return None

    # Get your RefSeq
    return soup_ul.li.text.split(' ')[1]


def refseq_to_assembly(refseq):
    '''Queries NCBI"s database to get the assembly code given a refseq.'''
    base_url = 'https://www.ncbi.nlm.nih.gov/assembly/REFSEQ/'

    if refseq is None:
        return None, None

    # Making a GET request
    sleep(1/2)
    get_request = requests.get(base_url.replace('REFSEQ', refseq))

    # Parse request
    soup = BeautifulSoup(get_request.content, 'html.parser')

    # Find object that contains relevant info
    soup_h1 = soup.find('h1', class_="marginb0 margin_t0")
    if soup_h1 is None:
        return None, None
    assembly = soup_h1.getText()

    # Extract real name of assembly's organism
    soup_h1 = soup.find(id="summary").find_all('dd')
    for elem_h1 in soup_h1:
        if "Organism name" in elem_h1.previous_sibling.get_text(strip=True):
            name = ' '.join(elem_h1.get_text(strip=True).split(' ')[:2])

    return assembly, name


def download_assembly_from_ncbi(refseq, assembly, type_, filename_out=None,
                                bar=True):

    '''Download the gunzipped sequence data from NCBI.'''

    def bar_custom(current, total):
        '''Produces a custom progress bar.'''
        print('\b'*999 +
              f"Downloadong: {current/total*100}%% [{current}/{total}] bytes",
              end='')

    def bar_none(*args, **kwargs):
        '''Empty progress bar'''
        return

    if filename_out is None:
        filename_out = f"{refseq}.fna.gz"

    if os.path.exists(filename_out):
        if VERBOSE:
            print(f'File {filename_out} already exists.')
            return True

    if assembly is None:
        assembly_name = refseq_to_assembly(refseq)
    else:
        assembly_name = assembly

    assembly_name = assembly_name.replace(' ', '_')
    header = refseq.split('_')[0]
    tag_1 = refseq.split('_')[1][0:3],
    tag_2 = refseq.split('_')[1][3:6]
    tag_3 = refseq.split('_')[1][6:9]

    post = '/'.join([header, tag_1, tag_2, tag_3, refseq])
    filename = f'{refseq}_{assembly_name}_{type_.lower()}_from_genomic.fna.gz'
    base_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
    url = f"{base_url}{post}_{assembly_name}/{filename}"

    if bar:
        bar = bar_custom
    else:
        bar = bar_none

    if os.path.exists(filename_out):
        os.remove(filename_out)

    try:
        wget.download(url, filename_out, bar=bar_custom)
        return filename_out
    except Exception as error:
        print('An exception occurred: {}'.format(error))
        return False


def sequence_from_assembly(query, refseq, assembly, type_):

    if (refseq is None) or (assembly is None):
        return '', ''

    filename_temp = f"{refseq}_data.temp"
    filename_assembly = f"./{type_}/{refseq}.fna.gz"

    # Create directory
    create_dir_if_doesnt_exist("./TYPE".replace("TYPE", type_))

    # Download assembly genomes, and look for the sequences of interest
    suc = download_assembly_from_ncbi(refseq,
                                      assembly,
                                      type_,
                                      filename_out=filename_assembly
                                      )
    if not suc:
        if VERBOSE:
            print(f'<W> Could not download the assembly code of {refseq}')
        return '', ''

    # Unzip file
    gunzip(filename_assembly, filename_temp)

    # Parse file searching for the query name
    header, body = get_sequence_from_fasta(filename_temp,  f"{query}")

    if os.path.exists(filename_temp):
        os.remove(filename_temp)

    return header, body


def get_sequence_from_fasta(filename: str, query: str):
    '''Reads a .fasta file and retrieves the sequences containing the query.'''

    header, body = '', ''

    with open(filename, 'r', encoding="utf-8") as file:
        line = 'X'
        while line != '':
            line = file.readline()

            if query in line:
                header = line.replace('\n', '')
                body = ''
                line = file.readline()
                while (len(line) > 0) and (line[0] != '>'):
                    body += line
                    line = file.readline()
                break
    return header, body


def format_header(header, query, name, truename, refseq, assembly, accession):
    '''Fasta sequence header standard formatting.'''
    loc = header_get_location(header)
    if truename is None:
        truename = name
    truename = truename.replace(' ', '_')
    header = f">{truename}:{query} [original_query={name}] " + \
             f"[refseq={refseq}] [assembly={assembly}] [seq_query={query}]" + \
             f" [loc={loc}] [accession={accession}]\n"

    return header


def header_get_location(header: str):
    '''Retrieves the location of a sequence from its header.'''

    if "location" in header:
        return header.split('location')[1].split("]")[0][1:]
    else:
        return "Unknown"


def gunzip(source_filepath, dest_filepath, block_size=65536):
    '''Unzips source_filepath into dest_filepath.'''

    import gzip

    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)


def create_dir_if_doesnt_exist(path):
    '''Creates a directory if it does not exist.'''

    if os.path.isfile(path):
        raise ValueError("Not a dir path.")

    if not os.path.isdir(path):
        os.mkdir(path)


def sequence_from_entrez_search(query):
    '''Returns a query sequence from a request to NCBIs entrez'''

    idlist = ncbi_esearch(query)
    fastas = ncbi_efetch(idlist[:5])

    if fastas is None:
        return '', ''

    for fasta in fastas:
        _h, _b = fasta
        if every_word_in_str(query.split(' '), _h):
            header = _h
            body = _b
            return header, body
    return '', ''


def every_word_in_str(query: list, reference: str):
    '''Checks if every word from a list is in a string.'''

    return np.all([word in reference for word in query])


def ncbi_esearch(term, database='nuccore', retmode='json'):
    '''Performs a search to NCBIs Esearch api'''

    url_search = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params_search = {
        'db': database,  # or nucleotide
        'term': term,
        'retmode': retmode,
        'sort': 'relevance'}

    sleep(1/2)
    response = requests.get(url_search, params=params_search, timeout=10)

    if response.status_code == 200:
        json = response.json()
        idlist = json['esearchresult']['idlist']
        return idlist
    else:
        print('<W> Error in the NCBI search.')
        return None


def ncbi_efetch(GI, db='nuccore', rettype='fasta', retmode='text'):
    '''Performs a query to NCBIs Efetch api'''

    if GI is None:
        return None

    url_efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params_fetch = {
        'db': db,  # or nucleotide
        'id': GI,
        'rettype': rettype,
        'retmode': retmode}

    sleep(1)
    response = requests.get(url_efetch, params=params_fetch, timeout=10)

    if response.status_code == 200:
        text = response.text.split('>')[1:]
        fastas = []
        for mixed in text:
            lines = mixed.split('\n')
            header = '>'+lines[0]
            body = '\n'.join(lines[1:])
            fastas.append([header, body])

        return fastas


def add_data(data, name=None, truename=None, refseq=None, assembly=None,
             accession=None, success=None):
    '''Adds a single entry to a dictionnary of lists.'''

    data["Name"].append(name)
    data["TrueName"].append(truename)
    data["Refseq"].append(refseq)
    data["Assembly"].append(assembly)
    data["Accession"].append(accession)
    data["SuccessWith"].append(None)
    return data


def is_sequence_in_fasta(filename, query_list):
    '''Checks if sequences is present in fasta file.'''

    if not os.path.isfile(filename):
        return False

    with open(filename, 'r', encoding="utf-8") as file:
        for line in file.readlines():
            if every_word_in_str(query_list, line):
                return True
    return False


def export_failed(organism, filename):
    '''Appends filename some string.'''

    with open(filename, 'a+', encoding="utf-8") as file:
        file.write(organism + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                prog='get_sequences_from_ncbi.py',
                description="get_sequences_from_ncbi.py belongs to the " +
                            "microguilds package. It searches for Refseqs " +
                            "and assembly codes in NCBI's database from " +
                            "organism names.",
                epilog="Please visit github.com/pyubero/microguilds for " +
                       "more details.")

    parser.add_argument('input', metavar='filename')
    parser.add_argument('-q', '--query',
                        metavar='sequence_name',
                        help="Name of genes or sequences queried.")
    parser.add_argument('-t', '--type',
                        metavar="type_",
                        help="Either CDS or RNA",
                        choices=["CDS", "RNA"])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help="Inline print out of detailed progress.")
    parser.add_argument('--clear_all',
                        action='store_true',
                        help='Clears all output before running the script.')

    # For interactive execution
    args = parser.parse_args(
        "organisms_input_clean.csv -q LepA -t CDS --verbose".split(' ')
        )

    # External variables
    FILENAME_IN = args.input
    QUERY = args.query
    TYPE = args.type
    CLEAR_ALL = args.clear_all
    VERBOSE = args.verbose

    # Internal variables
    OUT_CORRECT = './data/organismlist_correct.csv'
    OUT_FAILED = './data/organismlist_failed.csv'
    OUT_FASTA = f'./data/{QUERY}_sequences.fasta'
    nerrors = 0
    data = {"Name": [],
            "TrueName": [],
            "Refseq": [],
            "Assembly": [],
            "Accession": [],
            "SuccessWith": []}

    # Decrlate verbosity
    if VERBOSE:
        def verboseprint(*args, **kwargs):
            '''Prints if VERBOSE is true'''
            print(*args, **kwargs)
    else:
        def verboseprint(*args, **kwargs):
            '''Empty function'''
            return

    if CLEAR_ALL:
        if os.path.exists(OUT_CORRECT):
            os.remove(OUT_CORRECT)
        if os.path.exists(OUT_FASTA):
            os.remove(OUT_FASTA)
        if os.path.exists(OUT_FAILED):
            os.remove(OUT_FAILED)

    # 1. Load organism names
    orgnames = pd.read_csv(FILENAME_IN, header=0)["Name"]
    verboseprint(f"Found {len(orgnames)} input organisms.")

    # ... load previous session
    if os.path.exists(OUT_CORRECT):
        df_prevsession = pd.read_csv(OUT_CORRECT, keep_default_na=False)

        for key, _ in data.items():
            data[key] = df_prevsession[key].to_list()

        verboseprint(f"Loaded info on {len(df_prevsession)} " +
                     "organisms from previous session.")

    # ... for every organism in the list
    for jj, organism in enumerate(orgnames):

        verboseprint("--", jj, "-- ", organism)
        header, body = "", ""
        truename, refseq, assembly, accession = "", "", "", ""

        # Resume if already computed
        if organism in data["Name"]:
            idx = np.argwhere(organism == np.array(data["Name"])).flatten()[0]
            truename, refseq, assembly, accession = \
                [data[key][idx] for key in data][1:-1]

            if is_sequence_in_fasta(OUT_FASTA,
                                    [organism, refseq, assembly, accession]):
                verboseprint("Already found in ", OUT_CORRECT)
                verboseprint("")
                continue

        # 2. Obtain Refseq from name
        if refseq == "":
            refseq = name_to_refseq(organism)
        verboseprint('Refseq: ', refseq)

        # 3. Obtain assembly from refseq
        if assembly == "":
            assembly, truename = refseq_to_assembly(refseq)
        verboseprint('Assembly: ', assembly)
        verboseprint('Tru name: ', truename)

        # 4a. Search sequence(s) in assembly genome
        if (refseq != "") and (assembly != ""):
            header, body = sequence_from_assembly(QUERY,
                                                  refseq,
                                                  assembly,
                                                  TYPE
                                                  )
        verboseprint("Seq length: ", len(body.replace('\n', '')))

        # 4b. If empty, search in NCBI's nuccore
        if body == "":
            header, body = sequence_from_entrez_search(f"{organism} {QUERY}")
            accession = header.split(' ')[0][1:]
            verboseprint("Using EUTILS...")
            verboseprint("Accession: ", accession)
            verboseprint("Seq length: ", len(body.replace('\n', '')))

        # 5a. If failed...
        if body == "":
            nerrors += 1
            verboseprint('')
            export_failed(organism, OUT_FAILED)
            continue

        # 5b. If successful...
        header = format_header(header,
                               QUERY,
                               organism,
                               truename,
                               refseq,
                               assembly,
                               accession)

        with open(OUT_FASTA, 'a+', encoding="utf-8") as f:
            f.write(header)
            f.write(body)

        # 6. Update output table
        data = add_data(
            data,
            name=organism,
            truename=truename,
            refseq=refseq,
            assembly=assembly,
            accession=accession,
            success=None)

        pd.DataFrame.from_dict(data).to_csv(OUT_CORRECT, index=False)
        verboseprint('')

verboseprint(f"Found {nerrors} errors.")
