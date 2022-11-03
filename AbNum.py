################################################
# AbNum.py
################################################ 

"""Methods for numbering Ab amino acids with various schemes"""

import os
import pandas as pd
import re
from anarci import number
                      #CDR1  #FW2   #CDR2  #FW3   #CDR3  #FW4
CDR_idx = { "k_H":   ["31",  "36",  "50",  "66",  "95",  "103",  "113"],
            "c_H":   ["26",  "",    "52",  "57",  "95",  "103",  "113"],
            "i_H":   ["26",  "34",  "51",  "58",  "93",  "103",  "113"],
            "hum_H": ["26",  "36",  "50",  "66",  "93",  "103",  "113"],

                        #CDR1   #FW2   #CDR2  #FW3   #CDR3  #FW4
            "k_L":   ["24",  "35",  "50",  "57",  "89",  "98",  "109"],
            "c_L":   ["24",  "35",  "50",  "57",  "89",  "98",  "109"],
            "i_L":   ["27",  "33",  "50",  "53",  "89",  "98",  "109"],
            "hum_L": ["24",  "35",  "50",  "57",  "89",  "98",  "109"]}

def getNumsKabat(aaseq):
    """Get Kabat numbering and Ab sequence from ANARCI"""

    numbering, chain_type = number(aaseq, scheme="kabat")
    data = []

    for i in numbering:
        nums = str(i[0][0])+ i[0][1]
        nums = nums.replace(" ", "")
        data.append([chain_type, str(nums), i[1]])

    kabatdf = pd.DataFrame(data, columns = ["Chain", "KabatNum", "AA"])

    return kabatdf

def getNumsIMGT(aaseq):
    """Get Kabat numbering and Ab sequence from ANARCI"""

    numbering, chain_type = number(aaseq, scheme="imgt")

    data = []

    for i in numbering:
        nums = str(i[0][0])+ i[0][1]
        nums = nums.replace(" ", "")
        data.append([chain_type, str(nums), i[1]])

    kabatdf = pd.DataFrame(data, columns = ["Chain", "KabatNum", "AA"])

    return kabatdf

def defCDRs(scheme, kabatdf):
    """Find CDR indices"""

    scheme = scheme + "_" + kabatdf["Chain"][0] 

    #chothia weirdness
    if scheme == 'c_H':
        a35 = False
        b35 = False

        for index, row in kabatdf.iterrows():
            if row['KabatNum'] == '35A':
                a35 = True
            if row['KabatNum'] == '35B':
                b35 = True
        
        if a35 and b35:
           CDR_idx['c_H'][1] = "35"
        if a35:
            CDR_idx['c_H'][1] = "34"
        else:
            CDR_idx['c_H'][1] = "33"
        
    CDRs_idx = []

    #finding cdr starts & stops
    for index, row in kabatdf.iterrows():
        kabat_num = row[kabatdf.columns[1]]
        if kabat_num in CDR_idx[scheme]:
            CDRs_idx.append(index)

    count=0

    if len(CDRs_idx) == 5:
        last_idx = int(CDR_idx[scheme][5])
        for index,row  in kabatdf.iterrows():
            try:
                thisnum = int(row['KabatNum'])
                
            except:
                thisnum=0

            if thisnum > last_idx:
                if count == 0:
                    CDRs_idx.append(index)
                count+=1

    return CDRs_idx
