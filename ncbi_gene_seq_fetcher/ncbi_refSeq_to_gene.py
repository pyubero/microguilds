# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 18:25:32 2022

@author: Pablo
"""

import requests
import pycountry
import numpy as np
import pandas as pd
from time import sleep
import xmltodict as xml2dict
from matplotlib import pyplot as plt
from tqdm import tqdm
import xml.etree.ElementTree as ET




def refseq2assembly( refseq):
    BASE_URL = 'https://www.ncbi.nlm.nih.gov/assembly/REFSEQ/'
    
    # Making a GET request
    r = requests.get( BASE_URL.replace('REFSEQ', refseq))
    
    # Parse request
    soup = BeautifulSoup(r.content, 'html.parser')
    
    # Find object that contains relevant info
    s = soup.find('h1', class_="marginb0 margin_t0")
    if (s is None):
        print('Could not find anything for %s' % orgname )
        return None
    
    # Get your RefSeq
    return s.getText()    







nlines=0















