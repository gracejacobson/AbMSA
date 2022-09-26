########################################################
# BokehMSA
########################################################

"""Builds Bokeh plot for interactive antibody multisequence alignments"""

import numpy as np
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO, SeqIO
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.glyphs import Text, Rect
from AbNum import *

def getAbMSA(infile, scheme):
    """Main function for building plot"""

    #get alignment
    alnfile = getAln(infile)
    aln = AlignIO.read(alnfile, 'fasta')

    #get Ab kabat numbering
    kabatdict, idxs_dict = getNumsCdrs(infile, scheme)
    idxs_dict = correctCDRs(idxs_dict, aln)

    aln.sort(reverse=True)

    p = view_alignment(aln, idxs_dict, kabatdict)
    return p

########################################################
# Alignment & CDR Helper Functions
########################################################

def getNumsCdrs(infile, scheme):
    "Gets kabat numbers for alignment"
    idxs_dict = {}
    kabatdict = {}

    records = SeqIO.parse(infile, 'fasta')
    outfile = infile.replace(".fasta", "")

    for record in records:
        seq = str(record.seq)
        kabatdf = getNums(seq, outfile)
        CDRs_idx = defCDRs(scheme, kabatdf)
        kabatdict[str(record.id)] = [i for i in kabatdf['KabatNum']]
        idxs_dict[str(record.id)] = CDRs_idx

    return kabatdict, idxs_dict

def correctCDRs(idxs_dict, aln):
    """Defines CDR indices in alignment sequences"""
    for seq in aln:

        idxs = []
        a = str(seq.seq)
        id = seq.id

        for i in re.finditer('-', a):
            idxs.append(i.start())
        for z in idxs:
            for i in range(0, len(idxs_dict[id])):
                newlist = idxs_dict[id]
                if z <= idxs_dict[id][i]:
                    newlist[i] =  newlist[i] +1
                else:
                    newlist[i] = newlist[i]

                idxs_dict[id] = newlist

    return idxs_dict

def getAln(infile):
    """Gets alignment using clustam Omega"""
    outfile = infile.replace('.fasta','.aln')
    clustalomega_cline = ClustalOmegaCommandline(infile=infile, outfile=outfile, force=True)
    clustalomega_cline()
    return outfile

########################################################
# Bokeh Helper Functions
########################################################

def get_colors(seqs, consensus_seq):
    """make colors for bases in sequence"""
    consensus_seq = [i for i in consensus_seq]
    colors = []
   
    for sequence in seqs:
        for idx in range(0, len(consensus_seq)):
            if consensus_seq[idx] == '-':
                colors.append('white')

            else:
                if sequence[idx] != consensus_seq[idx]:
                    colors.append('white')
                else:
                    colors.append('blue')

    return colors

def getKabat(kabatdict, aln):
    """Array of kabat number for hover tool"""
    kabatnums = []

    for a in aln:
        count=0
        for idx in range(0, len(a.seq)):
            if a.seq[idx] == '-':
                kabatnums.append('N/A')
                count +=1
            else:
                kabatnums.append(kabatdict[a.id][idx-count])
    return kabatnums

def getCDR(aln, idxs_dict):
    """Array of Ab regions for hover tool"""
    rectx = []
    recty = []
    width = []
    y =0.5

    for rec in aln:
        a = str(rec.seq)
        id = str(rec.id)
        width1=idxs_dict[id][1]-idxs_dict[id][0]
        width2=idxs_dict[id][3]-idxs_dict[id][2]
        width3=idxs_dict[id][5]-idxs_dict[id][4]

        cdr1start = idxs_dict[id][0]+.5
        cdr2start= idxs_dict[id][2] +.5
        cdr3start= idxs_dict[id][4] +.5

        rectx.append(cdr1start + (width1/2))
        rectx.append(cdr2start + (width2/2))
        rectx.append(cdr3start + (width3/2))

        recty.append(y)
        recty.append(y)
        recty.append(y)
            
        width.append(width1)
        width.append(width2)
        width.append(width3)
        y+=1

    return rectx, recty, width

def getRegion(idxs_dict, aln):
    """Array of Ab regions for hover tool"""
    region = []
    for rec in aln:
        a = list(str(rec.seq))
        id = str(rec.id)
        if idxs_dict[id]:
            for i in range(0, len(a)):

                if i< idxs_dict[id][0]:
                    region.append("FWR1")
                if i>= idxs_dict[id][0] and i<idxs_dict[id][1]:
                    region.append("CDR1")
                if i>= idxs_dict[id][1] and i<idxs_dict[id][2]:
                    region.append("FWR2")
                if i>= idxs_dict[id][2] and i<idxs_dict[id][3]:
                    region.append("CDR2")
                if i>= idxs_dict[id][3] and i<idxs_dict[id][4]:
                    region.append("FWR3")
                if i>= idxs_dict[id][4] and i<idxs_dict[id][5]:
                    region.append("CDR3")
                if i>= idxs_dict[id][5]:
                    region.append('FWR4')
    return region

########################################################
# Bokeh Sequence Alignment View
########################################################

def view_alignment(aln,  idxs_dict, kabatdict, fontsize="9pt", plot_width=1300):
    """Builds MSA in Bokeh Plots"""

    #cutoffs for showing conservation
    summary = AlignInfo.SummaryInfo(aln)
    consensus_seq = summary.dumb_consensus(0.51, "-")

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    N = len(seqs[0])
    S = len(seqs)    

    #arrays for tool tips
    kabatnums = getKabat(kabatdict, aln) #kabat number for each AA
    colors = get_colors(seqs, consensus_seq) #color with >50% conservation
    region = getRegion(idxs_dict, aln) # CDR or FWR region
    aminoacid =  [z for rec in aln for z in list(range(0,len(rec.seq)))] # AA
    name = []
    for rec in aln:
        for z in range(0,len(rec.seq)):
            name.append(rec.id)

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5

    #array for CDR coords
    cdrx, cdry, cdrw = getCDR(aln, idxs_dict)

    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors, 
        aminoacid=aminoacid, kabatnums = kabatnums, name=name, region = region))
    plot_height = len(seqs)*15+125

    #view_range is for the close up view
    view_range = (0,130)
    tools="xpan, xwheel_zoom, save, reset"
    TOOLTIPS = [("ID","@name"),("AA", "@text"), ("Kabat Number","@kabatnums"), ('Region', '@region')]

    #sequence text view with ability to scroll along x axis
    p1 = figure(plot_width=plot_width, plot_height=plot_height,
                x_range=view_range, y_range=ids, tools=tools,
                min_border=0, toolbar_location='below')        

    source2= ColumnDataSource(dict(cdrx=cdrx, cdry=cdry, cdrw=cdrw))
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    rectCDR = Rect(x="cdrx", y="cdry", width="cdrw", height=1, fill_color=None, line_color ='red',
            fill_alpha=0, line_width=1.5)
    
    p1.add_glyph(source, glyph)
    r = p1.add_glyph(source, rects)
    p1.add_tools(HoverTool(tooltips=TOOLTIPS, renderers=[r]))
    p1.add_glyph(source2, rectCDR)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p1.title.text_font_size = "25px"

    return p1
