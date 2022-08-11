# myapp.py
import os
import io
import string
import tqdist
import pandas as pd
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool
import subprocess
import uuid
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import holoviews as hv
#import hvplot.pandas
import bokeh
import itertools
#import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import csv
import ete3
from ete3 import Tree
import PIL
from PIL import Image
from holoviews import opts
hv.extension('bokeh')
from bokeh.plotting import figure, show, output_file
from bokeh.models import (BasicTicker, BoxSelectTool, ColorBar, ColumnDataSource, TapTool, LinearColorMapper, PrintfTickFormatter, LinearInterpolator, Grid, ImageURL, LinearAxis, Plot, Range1d)
from bokeh.transform import transform
from bokeh.palettes import BrBG, PiYG, RdGy, RdYlGn, YlGnBu, Blues256, Viridis256
from random import random
import os, io, string, tqdist
import holoviews
from bokeh.layouts import row, column
from bokeh.models import Button, Slider, MultiSelect, Select, Div
from bokeh.models.widgets import FileInput
from bokeh.plotting import figure, curdoc
############################################# FUNCTIONS/VALS ######################################
#Specify path to step_04 final filtered exon directory and step_05 exon phylogeny directory
exon_path = './data/03c_final_exons/'
trees_path = './data/04_exon_phylogeny/'

#Calculate number of final filtered exons
exons = len(os.listdir(exon_path))

#Get a list of all taxa names that appear in final exon alignments
taxa_dict = list()
for fn in os.listdir(exon_path):
    path = os.path.join(exon_path, fn)
    for seq_record in SeqIO.parse(path, "fasta"):
        taxa_dict.append(seq_record.id)
#Filter taxon name list to keep only unique names
unique_names = sorted(list(set(taxa_dict)))

orth = []
for n in os.listdir(exon_path):
   orth.append(n)

############################################# PAGE LAYOUT SETUP ####################################
# add heading
heading = Div(text="""<font size="5"><b>Exon selection and visualization</b></font>""")


mds={'Dim1':[0.05,0.65], 'Dim2':[0.05,0.65], 'boot':[0, 0], 'taxa': [1, 1], 'size': [0.01, 0.01], 'exon': ["a", "b"], 'len': [0, 0]}
palette = list(reversed(Blues256))
source2 = ColumnDataSource(data=mds)
p2=figure(width=500, height=500, title="MDS of exon trees", toolbar_location="right", x_axis_location="below", tools = "tap, lasso_select, wheel_zoom, save, box_zoom,reset, hover, box_select", tooltips=[('exon', '@exon'), ('# of taxa', '@taxa'), ('bootstrap', '@boot')])
mapper2= LinearColorMapper(palette=palette, low=0, high=100)
p2.circle('Dim1', 'Dim2', source=source2, fill_alpha=0.8, radius='size',color=transform('boot', mapper2))
p2.xaxis.axis_label = 'Dim1'
p2.yaxis.axis_label = 'Dim2'
color_bar2 = ColorBar(color_mapper=mapper2, title="Bootstrap support")
p2.add_layout(color_bar2, 'right')
p2.select(TapTool).mode = 'append'

# add multiselect to select taxa to display
multi_select = MultiSelect(title="Select taxa to show:", options=unique_names, value=unique_names, height=480)
#select_og = Select(title="Select one outgroup:", options=unique_names, value=unique_names[0])
og_select = MultiSelect(title="Select outgroup:", options=unique_names, value=[unique_names[0]], size=4)
## add a button widget and configure with the call back
button = Button(label="Show plots", width=100, height=50, align="end")
heat_button=Button(label="Export exons", width=100, height=50, align="end")
mds_button=Button(label="Export exons from MDS plot", width=100, height=50, align="center")
root_button = Button(label="Root tree", width=100, height=50, align="end")
#Function to compute pairwise sequence similarity
def compute_similarity(seq_1, seq_2):
    #Make sure they are the same length
    if len(seq_1) != len(seq_2):
        raise ValueError('Sequences must be aligned and have the same length')
    seq_1 = seq_1.lower()
    seq_2 = seq_2.lower()    
    #Set up counters of length and similarity
    comp_length = 0
    num_sim = 0    
    #Iterate through each position in the sequences
    for base in range(len(seq_1)):        
        #Skip gaps
        if (seq_1[base] != '-') and (seq_2[base] != '-'):        
            #Increase the counter for compared length
            comp_length += 1            
            #Compare the two positions
            if seq_1[base] == seq_2[base]:                
                #Increase the similarity counter
                num_sim += 1                
    #Compute and return the percent similarity
    if comp_length == 0:
        #print('one (or both) of sequences is empty')
        score = 0
    elif num_sim == 0:
        score = 0   
    else:
        score = num_sim / comp_length        
    return score
    

#Pairwise seq similarity for all taxa; generates sp vs sp plot
def sim_by_sp(name):
    fn = os.path.join(exon_path, name)
    info=[]
    #length=[]
    with open(fn, 'r') as seq:
        ortho_file = SeqIO.parse(seq, 'fasta')
        #ortho_seqs = [record for record in ortho_file]
        ortho_seqs=[]
        for record in ortho_file:
            length=len(record)
            if record.id not in multi_select.value:
                continue
            ortho_seqs.append(record)
            
        if len(ortho_seqs) > 2:
            miss=1-len(ortho_seqs)/len(multi_select.value)
            b=[]
            for pair in itertools.product(ortho_seqs, repeat=2):
                s = compute_similarity(*pair)
                b.append(s)
                mean = sum(b) / len(b)
            data = {'Exon': name, 'Similarity' : mean, 'Missingness' : miss, 'Length': length}
        else:
            miss=1-len(ortho_seqs)/len(multi_select.value)
            mean=0
            data = {'Exon': name, 'Similarity' : mean, 'Missingness' : miss, 'Length': length}
        info1=pd.DataFrame(data, index=[name])
    return info1
p = Pool()
#  start = time.time()
async_result = p.map(sim_by_sp, orth)
p.close()
p.join()
    
# concatenate results into a single pd.Series
results = pd.concat(async_result)
results['miss_bins'] = pd.qcut(results['Missingness'], q=10, precision=3, duplicates='drop')
results['sim_bins'] = pd.qcut(results['Similarity'], q=10, precision = 3, duplicates='drop')
raw_source=ColumnDataSource(results)
print(results.head(n=5))
ex_df = {'exon':results['Exon'].astype('str'), 'sim_bins':results['sim_bins'].astype('str'),'miss_bins': results['miss_bins'].astype('str'), 'len':results['Length']}
ex_source=ColumnDataSource(ex_df)

df3 = results.groupby(['miss_bins', 'sim_bins']).size()
df3a = df3.reset_index(level=[0,1])
df3a = df3a.rename(columns={df3a.columns[2]: 'count'})
df3a['miss_bins']=df3a['miss_bins'].astype('str')
df3a['sim_bins']=df3a['sim_bins'].astype('str')
df3a['count']=pd.to_numeric(df3a['count'])
datax={'sim_bins':df3a['sim_bins'].astype('str'),'miss_bins': df3a['miss_bins'].astype('str'), 'count': pd.to_numeric(df3a['count'])}
source1 = ColumnDataSource(datax)

p1=figure(width=500, height=500, title="Number of exons by missingness and similarity", x_range=list(np.unique(source1.data['miss_bins'])), y_range=list(np.unique(source1.data['sim_bins'])), toolbar_location="right", x_axis_location="below", tools = "wheel_zoom, reset, hover, save, tap, box_select", tooltips=[('# of exons', '@count'), ('similarity', '@sim_bins'), ('missing', '@miss_bins')])
mapper = LinearColorMapper(palette=list(Viridis256), low=source1.data['count'].min(), high=source1.data['count'].max())
p1.rect(x="miss_bins", y="sim_bins", width=1, height=1, source=source1, line_color=None, fill_color=transform('count', mapper))
color_bar = ColorBar(color_mapper=mapper, ticker=BasicTicker(desired_num_ticks=25), title="Number of exons")
p1.add_layout(color_bar, 'right')
p1.axis.axis_line_color = None
p1.axis.major_tick_line_color = None
p1.axis.major_label_text_font_size = "14px"
p1.axis.major_label_standoff = 0
p1.xaxis.major_label_orientation = 1.0
p1.xaxis.axis_label = 'Avg missingness'
p1.yaxis.axis_label = 'Avg similarity'

p1.select(BoxSelectTool).select_every_mousemove = False
p1.select(TapTool).mode = 'append'

###################### TREE PANEL ###################################################
p3 = Div(text="", width=400, height=575)
#####################################################################################

#Calculate pairwise quartet distances
def qrt_dist(tr):
    score = []
    for key in sub_trees:
        t1 = Tree(tr)
        t2 = Tree(sub_trees[key])
        common_leaves = set(t1.get_leaf_names()) & set(t2.get_leaf_names())
        t1.prune(common_leaves)
        t2.prune(common_leaves)
        a = t1.write()
        b = t2.write()
        qscore = tqdist.quartet_distance(a, b)
        score.append(qscore)
    return score



def update_data():
    orth = []
    for n in os.listdir(exon_path):
        orth.append(n)
    
    p = Pool()
    #  start = time.time()
    async_result = p.map(sim_by_sp, orth)
    p.close()
    p.join()
        
    # concatenate results into a single pd.Series
    results = pd.concat(async_result)
    if (results['Missingness'] == 0).all(): 
        results['miss_bins'] = '[0; 0]'
    else:
        results['miss_bins'] = pd.qcut(results['Missingness'], q=10, precision=3, duplicates='drop')
    results['sim_bins'] = pd.qcut(results['Similarity'], q=10, precision = 3, duplicates='drop')
    raw_source.data=results
    ex_df2 = {'exon':results['Exon'].astype('str'), 'sim_bins':results['sim_bins'].astype('str'),'miss_bins': results['miss_bins'].astype('str'), 'len':results['Length']}
    ex_source.data=ex_df2

    df3 = results.groupby(['miss_bins', 'sim_bins']).size()
    df3a = df3.reset_index(level=[0,1])
    df3a = df3a.rename(columns={df3a.columns[2]: 'count'})
    df3a['miss_bins']=df3a['miss_bins'].astype('str')
    df3a['sim_bins']=df3a['sim_bins'].astype('str')
    df3a['count']=pd.to_numeric(df3a['count'])
    print(df3a)
    datay={'sim_bins':df3a['sim_bins'].astype('str'),'miss_bins': df3a['miss_bins'].astype('str'), 'count': pd.to_numeric(df3a['count'])}
    source1.data = datay
    #mapper = LinearColorMapper(palette='Viridis256', low=source1.data['count'].min(), high=source1.data['count'].max())
    color_bar = ColorBar(color_mapper=mapper, ticker=BasicTicker(desired_num_ticks=12), title="Number of exons")
    p1.x_range.factors=list(np.unique(source1.data['miss_bins']))
    p1.y_range.factors=list(np.unique(source1.data['sim_bins']))
    p1.rect(x="miss_bins", y="sim_bins", width=1, height=1, source=source1, line_color=None, fill_color=transform('count', mapper))
    mapper.low = min(source1.data['count'])
    mapper.high = max(source1.data['count'])

    

def selection_change(attrname, old, new):
    ex=pd.DataFrame(data=ex_source.data)
    taxa = multi_select.value
    path = trees_path
    global sub_trees
    a=source1.data['miss_bins'][new]
    b=source1.data['sim_bins'][new]
    a1=a.values
    b1=b.values
    miss=[]
    for mr in range(len(a1)):
        miss_range=ex[(ex['miss_bins'].astype('str') == a1[mr]) & (ex['sim_bins'].astype('str') == b1[mr])].index.tolist()
        miss.append(miss_range)
    miss_flat = [item for sublist in miss for item in sublist]
    fin_list = list(set(miss_flat))
    fina_list=[i.split('.', 1)[0] for i in fin_list]
    final_list=[path + str(x) for x in fina_list]
    print(final_list)
    
    #find all raxml best trees and save them as strings to a list
    trees={}
    for dirs in os.listdir(path):
        ortho = os.path.join(path,dirs)
        for file in os.listdir(ortho):
            if ".support" in file:
                tree = os.path.join(ortho, file)
                with open(tree,"r") as f:
                    a = f.read()
                    trees[ortho] = a
            else:
                continue
    print('Total number of exon trees: ', len(trees))
    #get subset of tree based on the final selected list of exons
    #sub_trees = {key: trees[key] for key in final_list}
    sub_trees = {}
    for key in final_list:
        if key in trees:
            sub_trees[key] = trees[key]
        else:
            continue
    #mds_source=ColumnDataSource(sub_trees)
    #get avg bootstrap values for all trees
    avg_bootstrap = []
    taxa = []
    for key in sub_trees:
        t = Tree(sub_trees[key])
        total_taxa = len(t)
        node_bootstrap = []
        for node in t.traverse():
            node_bootstrap.append(node.support)
        total_bootstrap = sum(node_bootstrap)
        number_nodes = len(node_bootstrap)
        mean = round(float(total_bootstrap)/float(number_nodes),2)
        avg_bootstrap.append(mean)
        taxa.append(total_taxa)

    p = Pool()
    qrt_result = p.map(qrt_dist, sub_trees.values())
    p.close()
    p.join()
    
    #results = pd.concat(async_result)
    df_q = pd.DataFrame(qrt_result, columns=sub_trees.keys())
    print(df_q)
    #MDS stuff; reduce number of dimensions to 2
    mds = MDS(random_state=0, dissimilarity='precomputed')
    X_transform = mds.fit_transform(df_q)
    
    df_mds = pd.DataFrame(data=X_transform)
    print(df_mds)
    df_mds['boot'] = avg_bootstrap
    df_mds['taxa'] = taxa
    df_mds['size']=df_mds['taxa']/900
    df_mds=df_mds.apply(pd.to_numeric)
    #fina_list=[i.split('.', 1)[0] for i in fin_list]
    df_mds['exon'] =[i.split('/')[3] for i in list(df_q.columns)]
    df_mds = df_mds.rename(columns={df_mds.columns[0]: 'Dim1'})
    df_mds = df_mds.rename(columns={df_mds.columns[1]: 'Dim2'})
    print(df_mds.head(n=5))
    source2.data = df_mds
    #mapper2= LinearColorMapper(palette=palette, low=min(source2.data['boot']), high=max(source2.data['boot']))
    mapper2.low = min(source2.data['boot'])
    mapper2.high = max(source2.data['boot'])
    #p2.x_range=list(source2.data['Dim1'])
    #p2.y_range=list(source2.data['Dim2'])

#############################################
#############################################
def mds_change(attrname, old, new):
    dir = './HiMAP2_viz/static'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))
    tr_udi = uuid.uuid4().hex.lower()[0:10]
    tr_img_path="./HiMAP2_viz/static/tree" + tr_udi + ".png"
    new=source2.selected.indices
    ex_trees=list(source2.data['exon'][new])
    ex_tr=['./data/04_exon_phylogeny/' + s for s in ex_trees]
    speciestree_dir = './HiMAP2_viz/static'
    path_trees = os.path.join(speciestree_dir,"exon_trees.tre")
    log_err_file = os.path.join(speciestree_dir,"sp_tree.log")
    astralpath = "/home/oksanav/software/Astral/astral.5.7.8.jar"
    max_mem = 8000
    picked_trees = {key: sub_trees[key] for key in ex_tr}    
    output_speciestree = os.path.join(speciestree_dir,"sp_tree.tre")
    java_max_mem =  "-Xmx" + str(max_mem) + "M"
    with open(path_trees, 'w') as output_newick:
        for tr in picked_trees.values():
            output_newick.write(tr)
    #Get species_tree
    with open(log_err_file, "wb") as err_log:
        cmd = ["java",java_max_mem,"-jar",astralpath,"-i",path_trees,"-o",output_speciestree] 
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err_log)
        out = p.communicate()
    
    t=Tree(output_speciestree)
    
    t.render(tr_img_path)
    text = "<img src = \'" + tr_img_path + "\'>"
    p3.text = text
    
#############################################
#############################################
def root_tr():
    stat_dir = './HiMAP2_viz/static'
    files_in_directory = os.listdir(stat_dir)
    filtered_files = [file for file in files_in_directory if file.endswith(".png")]
    for file in filtered_files:
        path_to_file = os.path.join(stat_dir, file)
        os.remove(path_to_file)
    sp_tr = Tree('./HiMAP2_viz/static/sp_tree.tre')
    tr_udi = uuid.uuid4().hex.lower()[0:10]
    tr_img_path="./HiMAP2_viz/static/tree" + tr_udi + ".png"
    
    if len(og_select.value) == 1:
        og=og_select.value[0]
    else:
        print(len(og_select.value), og_select.value[:])
        tx = og_select.value[:]
        og = sp_tr.get_common_ancestor(tx)
    sp_tr.set_outgroup(og)
    sp_tr.render(tr_img_path, h=480, units="px")
    text = "<img src = \'" + tr_img_path + "\'>"
    p3.text = text

#############################################    
#############################################    
def export_exons():
    new=source1.selected.indices
    ex=pd.DataFrame(data=ex_source.data)
    raw=pd.DataFrame(data=raw_source.data)
    taxa=multi_select.value
    os.makedirs("./data/selected_exons", mode=0o777, exist_ok=True)
    out_path="./data/selected_exons"
    for f in os.listdir(out_path):
        os.remove(os.path.join(out_path, f))
    a=source1.data['miss_bins'][new]
    b=source1.data['sim_bins'][new]
    a1=a.values
    print(a1)
    b1=b.values
    print(b1)
    miss=[]
    sim=[]
    final_list=[]
    for mr in range(len(a1)):
        miss_range=ex[(ex['miss_bins'].astype('str') == a1[mr]) & (ex['sim_bins'].astype('str') == b1[mr])].index.tolist()
        miss.append(miss_range)
    miss_flat = [item for sublist in miss for item in sublist]
    #sim_flat = [item for sublist in sim for item in sublist]
    final_list = set(miss_flat)
    print(len(miss_flat))
    #print(len(sim_flat))
    print(final_list)
    
    sub_raw=raw[raw["Exon"].isin(miss_flat)].drop('Exon', 1)
    sub_raw.to_csv('./data/selected_exons_info.csv', index = False)
    
    for ex in final_list:
        seqs=[]
        path = os.path.join(exon_path, ex)
        for seq_record in SeqIO.parse(path, "fasta"):
            if seq_record.id in taxa:
                seqs.append(seq_record)
            else:
                continue
        out_file = os.path.join(out_path, ex)
        SeqIO.write(seqs, out_file, "fasta")
    
#############################################    
#############################################    
def export_mds_exons():
    new=source2.selected.indices
    ex=pd.DataFrame(data=source2.data)
    raw=pd.DataFrame(data=raw_source.data)
    taxa=multi_select.value
    os.makedirs("./data/selected_mds_exons", exist_ok=True)
    ex_path=exon_path[:-1]
    out_path="./data/selected_mds_exons"
    for f in os.listdir(out_path):
        os.remove(os.path.join(out_path, f))
    a=source2.data['exon'][new]
    print(a)
    print(ex_path)
    a1=a
    
    fin_list=list(a1)
    #sim_flat = [item for sublist in sim for item in sublist]
    final_list = set(fin_list)
    print(final_list)
    fina_list=[s + '.fasta' for s in fin_list]
    sub_raw=raw[raw["Exon"].isin(fina_list)].drop('Exon', 1)
    sub_raw.to_csv('./data/selected_mds_exons_info.csv', index = False)
    
    for ex in final_list:
        seqs=[]
        #ex1=ex.split('/')[3]
        #print(ex1)
        path = os.path.join(ex_path, ex + ".fasta")
        for seq_record in SeqIO.parse(path, "fasta"):
            if seq_record.id in taxa:
                seqs.append(seq_record)
            else:
                continue
        out_file = os.path.join(out_path, ex + ".fasta")
        SeqIO.write(seqs, out_file, "fasta")


#update heatmap
button.on_click(update_data)
#heatmap selection action
source1.selected.on_change('indices', selection_change)
#tree view update when mds tree points are selected
source2.selected.on_change('indices', mds_change)
#heatmap exon export
heat_button.on_click(export_exons)
#mds exon export 
mds_button.on_click(export_mds_exons)
#re-root tree
root_button.on_click(root_tr)

# set layout and add to the document
curdoc().add_root(column(heading, row(column(multi_select, button), column(p1, heat_button), column(p2, mds_button), column(row(p3), row(og_select, root_button)))))
