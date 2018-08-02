
# coding: utf-8

# # Modeling of gene expression dynamics

# ### Imports

# In[10]:

print('importing modules and custom classes')
# math
import numpy as np
np.random.seed(20170601)
# statistics
import scipy as sp
# data management
import pandas as pd
# plotting
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
# network analysis
import networkx as nx
# utility
from copy import deepcopy
# time
import time

# evo devo
from evodevo import *
from evodevo import Organism as org
print('imports complete')


# ### Function to generate and characterize a gene regulatory network in an 'organism'

# In[12]:

# function to generate an organism with a set of characters
def generate_organism(n, comp, a, mxi, epsilon):    
    print('generating a stable gene regulatory network\nGRN paramters:\nn=',n,
          'comp=',comp,'a=',a,'mxi=',mxi,'epsilon=',epsilon)
    i = 1
    while True:
        
        timestamp = time.time()
        name = str('fndr.' + str(i) + '.' + str(timestamp))
        founder = org(name)
        founder.wc = 10
        founder.make_grn(n, comp, a, mxi)
        founder.start_development(founder.s0, epsilon)

        if not founder.stabilized:
            del name, founder, timestamp
            i+=1
        else:            
            # how large of the space to sample
            chunk_size = 2**founder.n
            gob_id = []
            gob_id = np.arange(0,1)
            
            # boolean options, yes, no for each gene
            options = [0.0, 1.0]
            
            # creates an index for each chunk
            gob_idx = []
            for gob in gob_id:
                gob_idx = gob_id[gob] * chunk_size

            # explore attractor space
            vectors, states, paths, attractors, sbuff  = generate_combinations(founder.n, options, gob_idx, chunk_size, founder)

            # filtering candidates based on path length and number of final fixed points
            # the goal is to have networks with more than 2 states, every state has a null state
            # this will catch any network that produces at least two active final fixed points
            if np.min(paths[1:]) >= founder.wc and np.max(paths[1:]) < mxi-1 and len(np.unique(states, axis=0)) > 2:
                try:
                    # something to try: make many cells, etc...
                    print('genesis complete')
                    paths  = np.subtract(paths,founder.wc-1)
                    print('min', np.min(paths[1:]), 'window', founder.wc, 'max', np.max(paths[1:]), 'mean', np.mean(paths), 'var',np.var(paths))
                    return founder, vectors, states, paths, attractors, sbuff
                except:
                    raise Exception('Failed to find a grn with stablizing states')
                    i+=1  
            else:
                i+=1
print('ready for parameters')


# ### Set parameters

# In[13]:

# set gene regulatory network parameters
n, comp, a, mxi, epsilon = 10, 0.5, 100, 100, 0.0001
print('grn paramters are:', n, comp, a, mxi, epsilon)


# ### Generate an organism and characterize the state space

# In[14]:

# make an organism
f2, fvectors, fstates, fpaths, fattractors, fsbuff = generate_organism(n,comp,a,mxi,epsilon)


# ### Find unique attractors and final output states (cell types)

# In[15]:

# processing state space information
a = np.unique(fattractors, axis=0)
b = np.unique(fstates, axis=0)
print('# of attractors:', len(a),'\n# of steady states', len(b))


# ### Characterization and description of Gene Regulatory Network Properties

# In[16]:

# sort attractor states and fixed points against the vector reference
fvectors = np.asarray(fvectors).squeeze()

print('finding indices of network state space')
# indices of unique fixed points and attractors states
fpidx = find_state_index(np.unique(fstates,axis=0),fvectors)
faidx = find_state_index(np.unique(fattractors,axis=0),fvectors)

# indices of for all states
all_fpidx = find_state_index(fstates,fvectors)
all_faidx = find_state_index(fattractors,fvectors)

# make simple tuples of fixed points and attractors
ufps = np.unique(all_fpidx)
ufps = tuple(x for x in ufps)
uatts = np.unique(all_faidx)
uatts = tuple(x for x in uatts)

# pull out the subnetworks
idxs = np.arange(0, 2**n)
full_net = [(idxs[x],all_faidx[x], all_fpidx[x]) for x in idxs] # full network space
att_net = [(idxs[x],all_faidx[x]) for x in idxs] # attractor network space
all_paths = [(idxs[x],fpaths[x]) for x in idxs] # all path lengths in tuple form
fp_net = [(idxs[x],all_fpidx[x]) for x in idxs] # fixed point steady state network space

print('index of unique fixed point states', fpidx)
print('index of unique attractor states', faidx)


# ### Plot W matrix

# In[17]:

# plotting network matrix
print('plotting gene regulatory network as a matrix')
# plot original founder properties as example:
fig, ax = plt.subplots(figsize=(5,4), ncols=1, dpi=300)
w=f2.w
cax = plt.imshow(w, origin='lower', interpolation='none', cmap=plt.cm.get_cmap('PRGn', 9), aspect='equal', vmin=-2, vmax=2)
plt.colorbar(cax, fraction=0.046, pad=0.04)

plt.title('Gene Regulatory Network')
plt.setp(ax,yticks=np.arange(0,n))
plt.setp(ax,xticks=np.arange(0,n))
plt.xlabel('binding affinity')
plt.ylabel('genes')
plt.tight_layout()
plt.savefig('w.png')
plt.show()
plt.close()


# ### Plot network view

# In[18]:

print('plotting classic network view with nodes and edges for gene regulatory network')
# function to draw network using the position and sorting activating vs repressing
def draw_network(G, pos, ax, sg=None):
    nx.draw_networkx_nodes(G,pos, nodelist=G.nodes, node_size=120, node_color='k',alpha=.9)
    nx.draw_networkx_nodes(G,pos,node_size=100, node_color='w',alpha=1)
    for n in G:
        c=Circle(pos[n],radius=.01,alpha=0.1)
        ax.add_patch(c)
        G.node[n]['patch']=c
        x,y=pos[n]

    nx.draw_networkx_labels(G, pos,font_size=10, font_family='sans-serif', font_color='k')

    seen={}
    for (u,v,d) in G.edges(data=True):
        n1=G.node[u]['patch']
        n2=G.node[v]['patch']
        rad=0.1
        if (u,v) in seen:
            rad=seen.get((u,v))
            rad=(rad+np.sign(rad)*0.1)*-1
        colors = ['green', 'purple', 'black', 'white']
        lw = 2

        if d['weight']>0.0:
            alpha=0.5       
            e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                                arrowstyle='-|>',
                                connectionstyle='arc3,rad=%s'%rad,
                                mutation_scale=10.0,
                                lw=lw,
                                alpha=alpha,
                                color=colors[0])
            seen[(u,v)]=rad
            ax.add_patch(e)
        else:
            alpha=0.25     
            e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                                arrowstyle='-|>',
                                connectionstyle='arc3,rad=%s'%rad,
                                mutation_scale=10.0,
                                lw=lw,
                                alpha=alpha,
                                color=colors[1])
            seen[(u,v)]=rad
            ax.add_patch(e)           
    return e 

fig, ax = plt.subplots(figsize=(3,3), ncols=1, dpi=300)
G = nx.from_numpy_array(np.asarray(w))
pos=nx.circular_layout(G) # positions for all nodes
G2 = draw_network(G, pos, ax)
plt.title('GRN (W) Network View')
plt.axis('off')
plt.savefig('w_network.png')
plt.show()
plt.close(fig)


# ### Plot state space

# In[20]:

print('plotting state space sets')
# plot state space matrices
fig2 = plt.figure(figsize=(8,11),dpi=600)
ax = fig2.add_subplot(211, aspect='auto')
cax = plt.imshow(fstates.T, origin='upper', interpolation='none', cmap=plt.cm.get_cmap('YlGnBu_r', 2), aspect='auto',vmin=0,vmax=1, alpha=0.8)
plt.colorbar(cax, ticks = [0,1])
plt.title('all terminal cell states')
plt.setp(ax,yticks=np.arange(0,n))

ax2 = fig2.add_subplot(212, aspect='auto')
cax = plt.imshow(np.transpose(fattractors), origin='upper', interpolation='none', cmap=plt.cm.get_cmap('YlGnBu_r', 2), aspect='auto',vmin=0,vmax=1, alpha=0.8)
plt.colorbar(cax, ticks = [0,1])
plt.title('attractor basin network space')
plt.setp(ax2, yticks=np.arange(0,n))
plt.tight_layout()
plt.savefig('state_space.png')
plt.show()
plt.close(fig2)    

fig3 = plt.figure(figsize=(8,3),dpi=600)

ax3 = fig3.add_subplot(121, aspect='auto')
cax = plt.imshow(np.unique(fstates,axis=0).T, origin='upper', interpolation='none', cmap=plt.cm.get_cmap('YlGnBu_r', 2), aspect='auto',vmin=0,vmax=1, alpha=0.8)
plt.colorbar(cax, ticks = [0,1])
plt.title('terminal cell states')
plt.setp(ax3,yticks=np.arange(0,n))    

ax4 = fig3.add_subplot(122, aspect='auto')
cax = plt.imshow(np.transpose(np.unique(fattractors,axis=0)), origin='upper', interpolation='none', cmap=plt.cm.get_cmap('YlGnBu_r', 2), aspect='auto',vmin=0,vmax=1, alpha=0.8)
plt.colorbar(cax, ticks = [0,1])

plt.title('unique attractor states')
plt.setp(ax4,yticks=np.arange(0,n))

plt.tight_layout()
plt.savefig('unique_states.png')
plt.show()
plt.close(fig3)    


# ### Plot path length distribution

# In[21]:

# use a histogram to plot the distribution of path lengths to steady state/stability
fig = plt.figure(figsize=(6,4), dpi=300)
ax = fig.add_subplot(311, aspect='auto')
plt.hist(fpaths*np.less(fpaths,mxi-2))
plt.xlabel('time to stability (path length)')
plt.ylabel('frequency')
plt.savefig('paths.png')
plt.show()
plt.close()


# ### Plot and analyze attractor network

# In[23]:

# analyze attractor network data
print('analyzing attractor network space and graphing states as nodes in a spring layout (kamada-kawai)')
# generate a graph using the attractor network computed
G = nx.Graph(att_net)

# create the path length cost function layout
pos = nx.kamada_kawai_layout(G)

# characterize the nodes of attractor space by degree and sort out states
degs = nx.degree(G)

ddf = pd.DataFrame.from_records(degs, columns=['node', 'degree'])
ddf = ddf.sort_values(by=['degree'])

# finding garden of eden by degree
garden_of_eden = list([x for x in ddf['degree'] if x ==1])
print('# garden of eden states', len(garden_of_eden))

# finding attractors by degree > 1
attractor_basin = list([x for x in ddf['degree'] if x >1])
print('# num of attractors', len(attractor_basin))

# finding strong attractors by degree > 25
strong_attractors = list([x for x in ddf['degree'] if x > 20])
print('# num of strong attractors', len(strong_attractors))

# finding fixed points and comparing with list of strongest attractors by degree
print('\ntop attractors by degree')
print(ddf['degree'][-10:])

# attractor network graphing
fig = plt.figure(figsize=(8,8), dpi=600)
ax = fig.add_subplot(111, aspect='auto')
ax.axis('off')
print('plotting')
# draw all states as small and black
nx.draw_networkx_nodes(G, pos, node_size=1.4, node_color='k', alpha=0.2)
# draw unique attractors points medium and blue
nx.draw_networkx_nodes(G, pos, node_size=5, nodelist=uatts, node_color='b', alpha=0.4)
# draw unique fixed points large and red
nx.draw_networkx_nodes(G, pos, node_size=7.5, nodelist=ufps, node_color='r', alpha=0.8)
# draw semitranparent network edges
nx.draw_networkx_edges(G, pos, alpha=0.2)
# draw small text node state labels
nx.draw_networkx_labels(G, pos, font_size=2)
plt.savefig('attractors.png')
plt.show()
plt.close()


# # Task 1: Alter the connectivity (c/comp) from .14 to .75 and plot the results.
# ##### hint: n, c, a, mxi, epsilon = 10, ?, 100, 100, 0.00001

# ____

# # Task 2: Rerun simulation using 7, 8, or 9 genes
# ##### hint: n, c, a, mxi, epsilon = ?, .5, 100, 100, 0.0001

# ____

# # Task 3: Plot the gene expression dynamics in the fsbuff variable
# ##### hint: reshape your fsbuff contents and use imshow to visualize, remember to save the graphics
# ##### warning do NOT save all 1024 gene expression sets, this will use all of your disk space and quota on pythonanywhere!!!

# In[ ]:



