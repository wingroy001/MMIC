from __future__ import division
import networkx as nx
import os,random
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import Counter
from itertools import chain
import math
from scipy.special import comb, perm
import time
from random import choice
p=0.1#边的权重
G1= nx.read_edgelist('D:\\soc-sign-bitcoinalpha1.txt', 
  create_using=nx.DiGraph())  #读取图1
G2= nx.read_edgelist('D:\\soc-sign-bitcoinalpha2.txt', 
  create_using=nx.DiGraph())  #读取图2
for u, v, d in G1.edges(data=True):  #G1边的权
    d['weight'] = p
for u, v, d in G2.edges(data=True):  #G2边的权
    d['weight'] = p    
#nx.write_edgelist(G1, 'D:\\111.txt',comments="#", delimiter='  ', encoding='utf-8')
#G3=nx.read_edgelist("D:\\11.txt",create_using=nx.DiGraph())
G1.add_nodes_from(G2.nodes(data=True))
nodelist = list(G1.nodes())
nodelist.sort()
mapping = {old_label:new_label for new_label, old_label in enumerate(nodelist)}
H1= nx.relabel_nodes(G1, mapping)

G2.add_nodes_from(G1.nodes(data=True))
nodelist = list(G2.nodes())
n1=len(G2.nodes())
nodelist.sort()
mapping = {old_label:new_label for new_label, old_label in enumerate(nodelist)}
H2= nx.relabel_nodes(G2, mapping)

mapping1={i:i+n1 for i in range(n1)}
H2= nx.relabel_nodes(H2, mapping1)
p1=0.1  #加跨层边的权
E1=[(i,i+n1,p1) for i in range(n1)]  #添加跨层边及其权
E2=[(i+n1,i,p1) for i in range(n1)]
E=E1+E2
G=nx.compose(H1,H2)
G.add_weighted_edges_from (E)
#print G.get_edge_data(4, 0)   #边（1,6）的权
#plt.figure(4) 
#nx.draw(G,with_labels=True,node_size=800)
nx.write_weighted_edgelist(G, 'D:\\soc-sign-bitcoinalpha0.txt')
print nx.info(G)
print "done"
#生成的图写入文档0.txt                  
#G0= nx.read_edgelist('D:\\0.txt', nodetype=int,data=(('weight',float),), create_using=nx.DiGraph())
#print nx.info(G)
#读取图的方式
#G = nx.read_edgelist('D:\\00.txt', nodetype=int, data=(('weight',float),),create_using=nx.DiGraph())
#plt.figure(4) 
#nx.draw(G,with_labels=True,node_size=8)
