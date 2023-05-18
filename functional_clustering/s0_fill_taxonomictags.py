import numpy as np
import pandas as pd
import requests
import warnings
from bs4 import BeautifulSoup
from time import sleep
from tqdm import tqdm


def name_to_eolcode(orgname):
    '''Searches EOL for an organisms code given its name.'''

    orgname = orgname.replace(' ', '+')
    url = f"https://eol.org/search?utf8=%E2%9C%93&q={orgname}"
    classname = "uk-link-reset uk-display-block"

    sleep(DELAY)
    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.Timeout:
        print(f">>> Timeout <<< Please check {url}")
        return name_to_eolcode(orgname)


    if response.status_code > 200:
        warnings.warn(f"Request status code {response.status_code}.")
        return None

    soup = BeautifulSoup(response.text, 'html.parser')
    elements = soup.find_all("a", {"class": classname})

    if len(elements) == 0:
        return None

    # elif len(elements) > 1:
    #     warnings.warn(
    #         f'''Found more than one entry for {orgname} in EOL db.
    #         Check {response.url} for more info.'''
    #     )

    return int(elements[0].get('href').split('/')[-1])


def eolcode_to_taxonomy(eolcode, boxname="NCBI"):
    '''Searches EOL for an organisms taxonomy given its code.'''

    if eolcode is None:
        return None

    url = f"https://eol.org/pages/{eolcode}/names"

    sleep(DELAY)
    try:
        response = requests.get(url, timeout=10)
    except requests.exceptions.Timeout:
        print(f">>> Timeout <<< Please check {url}")
        return eolcode_to_taxonomy(eolcode, boxname=boxname)

    if response.status_code > 200:
        warnings.warn(f"Request status code {response.status_code}.")
        return None

    soup = BeautifulSoup(response.text, 'html.parser')

    # Find bacterium name
    orgname = soup.find_all("div", {"class": "names"})[0].h1.text

    # Find bacterium taxonomy
    elements = soup.find_all("h4", {"class": "ui header"})
    element = [elem for elem in elements if boxname in elem.text]
    if len(element) == 0:
        warnings.warn(f"Could not find {boxname} taxonomy of {orgname}.")
        return None

    values = []
    for elem in element[0].find_next_sibling("div").find_all("a")[2:]:

        if elem.text in orgname:
            break
        if "group" in elem.text:
            continue
        else:
            values.append(elem.text)
    return values






taxonomy = values


eolcode = 11683786


FILENAME = "./data/organismlist_correct.csv"
DELAY = 0.8
TAGS = {
    "Domain": "d_",
    "Phylum": "p_",
    "Class": "c_",
    "Order": "o_",
    "Family": "f_",
    "Genus": "g_",
    "Species": "s_"
}

data = pd.read_csv(FILENAME)
data = data.fillna(np.nan)

TAXONOMIC_COLUMNS = []
for ii, row in tqdm(data.iterrows(), total=len(data)):
    name = row["Name"]
    code = name_to_eolcode(name)
    taxonomy = eolcode_to_taxonomy(code, boxname="NCBI")

    if taxonomy is None:
        taxonomy = eolcode_to_taxonomy(code, boxname="EOL Dynamic Hierarchy")

    if taxonomy is not None:
        if len(taxonomy) < 5:
            warnings.warn(f"Uncomplete taxonomy for {name}, {code}.")
            taxonomy = [*taxonomy, None, None, None, None, None][:5]

        TAXONOMIC_COLUMNS.append([*taxonomy, *name.split(' ')])

    else:
        TAXONOMIC_COLUMNS.append([None,]*7)
        warnings.warn(f"Did not find anything for {name}, {code}.")



for jj, key in enumerate(TAGS):
    data[key] = [tx[jj] for tx in TAXONOMIC_COLUMNS]


data.to_csv("test.tsv", sep='\t', index=False)
