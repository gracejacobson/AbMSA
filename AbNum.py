################################################
# AbNum.py
################################################ 

"""Methods for numbering Ab amino acids with various schemes"""

import pandas as pd
from anarci import number

                      #CDR1  #FW2   #CDR2  #FW3   #CDR3  #FW4
CDR_idx_kabat = { "k_H":   ["31",  "36",  "50",  "66",  "95",  "103",  "113"],
                  "c_H":   ["26",  "",    "52",  "57",  "95",  "103",  "113"],
                  "i_H":   ["26",  "34",  "51",  "58",  "93",  "103",  "113"],

                            #CDR1   #FW2   #CDR2  #FW3   #CDR3  #FW4
                  "k_L":   ["24",  "35",  "50",  "57",  "89",  "98",  "109"],
                  "c_L":   ["24",  "35",  "50",  "57",  "89",  "98",  "109"],
                  "i_L":   ["27",  "33",  "50",  "53",  "89",  "98",  "109"]}

                           #CDR1  #FW2   #CDR2  #FW3   #CDR3   #FW4
CDR_idx_IMGT = { "k_H":   ["36",  "41",  "55",  "75",  "107",  "118",  "128"],
                 "c_H":   ["27",  "38",  "57",  "65",  "107",  "118",  "128"],
                 "i_H":   ["27",  "39",  "56",  "66",  "105",  "118",  "128"],

                            #CDR1   #FW2   #CDR2  #FW3   #CDR3  #FW4
                 "k_L":   ["24",  "41",  "56",  "70",  "105",  "118",  "128"],
                 "c_L":   ["24",  "41",  "56",  "70",  "105",  "118",  "128"],
                 "i_L":   ["27",  "39",  "56",  "66",  "105",  "118",  "128"]}

def getNumsKabat(aaseq):
    """Get Kabat numbering and Ab sequence from ANARCI"""

    numbering, chain_type = number(aaseq, scheme="kabat")
    data = []

    for i in numbering:
        nums = str(i[0][0])+ i[0][1]
        nums = nums.replace(" ", "")
        data.append([chain_type, str(nums), i[1]])

    df = pd.DataFrame(data, columns = ["Chain", "Num", "AA"])
    return df

def getNumsIMGT(aaseq):
    """Get Kabat numbering and Ab sequence from ANARCI"""

    numbering, chain_type = number(aaseq, scheme="imgt")
    data = []

    for i in numbering:
        nums = str(i[0][0])+ i[0][1]
        nums = nums.replace(" ", "")
        data.append([chain_type, str(nums), i[1]])

    df = pd.DataFrame(data, columns = ["Chain", "Num", "AA"])
    return df

def defCDRs(scheme, df, numbering):
    """Find CDR indices"""

    scheme = scheme + "_" + df["Chain"][0] 

    if numbering == "k":
        #chothia weirdness
        if scheme == 'c_H':
            a35 = False
            b35 = False

            for index, row in df.iterrows():
                if row['Num'] == '35A':
                    a35 = True
                if row['Num'] == '35B':
                    b35 = True
            
            if a35 and b35:
                CDR_idx_kabat['c_H'][1] = "35"
            if a35:
                CDR_idx_kabat['c_H'][1] = "34"
            else:
                CDR_idx_kabat['c_H'][1] = "33"
            
        CDRs_idx = []

        #finding cdr starts & stops
        for index, row in df.iterrows():
            kabat_num = row[df.columns[1]]
            if kabat_num in CDR_idx_kabat[scheme]:
                CDRs_idx.append(index)
    
    if numbering == "i":
        CDRs_idx = []

        #finding cdr starts & stops
        for index, row in df.iterrows():
            kabat_num = row[df.columns[1]]
            if kabat_num in CDR_idx_IMGT[scheme]:
                CDRs_idx.append(index)

    return CDRs_idx
