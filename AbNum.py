################################################
# AbNum.py
################################################ 

"""Methods for numbering Ab amino acids with various schemes"""

import os
import pandas as pd
import re
                          #CDR1  #FW2   #CDR2  #FW3   #CDR3  #FW4
CDR_idx = { "k_H":   ["31",  "36",  "50",  "66",  "95",  "103",  "113"],
            "c_H":   ["26",  "",    "52",  "57",  "96",  "102",  "113"],
            "m_H":   ["30",  "36",  "47",  "59",  "93",  "102",  "113"],
            "i_H":   ["26",  "36",  "50",  "57",  "93",  "103",  "113"],

                        #CDR1   #FW2   #CDR2  #FW3   #CDR3  #FW4
            "k_LK":   ["24",  "35",  "50",  "57",  "89",  "98",  "109"],
            "c_LK":   ["26",  "33",  "50",  "53",  "91",  "97",  "109"],
            "m_LK":   ["30",  "37",  "47",  "56",  "89",  "97",  "109"],
            "i_LK":   ["27",  "33",  "50",  "52",  "89",  "98",  "109"]}

def getNums(aaseq, outfile):
    """Get Kabat numbering and Ab sequence from ANARCI"""

    #calling anarci
    os.system('ANARCI -i ' + aaseq + ' -s k -o ' + outfile + '.txt')
    
    #formatting raw df
    rawdf = pd.read_csv(outfile + '.txt', comment='#', header=None)
    rawdf.drop(rawdf.tail(1).index, inplace=True)

    data =[]
    for index, row in rawdf.iterrows():

        # chain, kabat num, and aa
        kabatList = re.findall(r"([HKL])\s(\d*)\s*([A-Z]?)\s*([A-Z-])",rawdf.at[index, rawdf.columns[0]] )
        
        if kabatList:
            if len(kabatList[0]) ==4:
                                #chain             #kabat num                              #amino acid
                data.append([kabatList[0][0], str(kabatList[0][1] + kabatList[0][2]), kabatList[0][3]])

            else:
                pass

    #formatting df
    kabatdf = pd.DataFrame(data, columns = ["Chain", "KabatNum", "AA"])

    return kabatdf

def defCDRs(scheme, kabatdf):
    """Find CDR indices"""

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

    return CDRs_idx