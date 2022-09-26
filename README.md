# AbMSA
This is a multi-sequence alignment viewer for antibody chains. Check out the Jupyter notebook for an example.

**Note: Requires ANARCI**. ANARCI is used to number amino acids (with Kabat) and get CDR definitions. 

This MSA viewer was based on code by [Damien Farrell](https://github.com/dmnfarrell) at UC Davis from [this](https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner) article. 

Each alignment requires a fasta file, and generates an interactive Bokeh figure. Along the way, an alignment file and .txt file is generated. A .txt file is overwritten with each chain.

Use the following schemes when inputting CDR definitions for getAbMSA(infile, scheme): 

Heavy Chain:
- "k_H": Kabat 
- "c_H": Chothia 
- "m_H": Martin 
- "i_H": IMGT 

Light Chain:
- "k_LK": Kabat
- "c_LK": Chothia
- "m_LK": Martin
- "i_LK": IMGT