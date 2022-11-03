# AbMSA

**REQUIRES CLUSTALO**

**Python libraries**: anarci, bokeh, biopython, panel, pandas.

This is a multi-sequence alignment viewer for antibody chains. Check out the Jupyter Notebook example using patented PD-L1 heavy chains!

![Example Alignment](/example/Example_Alignment.png)

***

This MSA viewer was based on code by [Damien Farrell](https://github.com/dmnfarrell) at UC Davis from [this](https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner) article. 

Each alignment requires a fasta file, and generates an interactive Bokeh figure. Along the way, an alignment file  is created using the BioPython clustalo wrapper.

Default Ab numbering is Kabat. You can specify what numbering scheme and what CDR definition to use: 
- "k": Kabat 
- "i": IMGT 
- "c": Chothia 
