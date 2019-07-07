from __future__ import division
import operator
import networkx as nx
import os,random
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import Counter
from itertools import chain
import math
from scipy.special import comb, perm
import time,copy
G= nx.read_edgelist('D:\\graph300.txt', nodetype=int,
  data=(('weight',float),), create_using=nx.DiGraph())
print '已经赋权'   
print nx.info(G) 
Q=G.reverse(copy=True)
n=len(G.nodes()) #图G的节点数n
m=len(G.edges()) #图G的边数m
ep,l=0.1,1 #种子数k，ep是参数epsilon,参数l 作为输入
acc=1-1/math.e-ep  #accuracy
pro=1-2/(n**l)#probability
print "accuracy is:",acc,"probability is:",pro  

m1=2  #网络的层数
n1=n/m1 #每层的点数 ，也就是实体的个数
def runIC (G, S): #运行一次IC模型
    T = deepcopy(S) # copy already selected nodes
    i = 0
    while i < len(T):
        for v in G[T[i]]: # for neighbors of a selected node
            if v not in T: # if it wasn't selected yet
                w = G[T[i]][v]['weight'] # count the number of edges between two nodes
                if random.random() <= w: # if at least one of edges propagate influence
                     T.append(v)
        i += 1
    return T
def F(R,S): #被S覆盖的RR集在总的RR集R中的比例
 b1=len(R)
 for j in range(len(S)):
  R=[R[i] for i in range(len(R)) if (S[j] not in R[i])] #没被S覆盖的RR集
 F=1-len(R)/b1
 return F
def NodeSelection(R,k):  #输出覆盖最多RR集的k个点
  S=set() #种子集S初始化
  B=R
  for j in range(k):  
    A= list(chain(*B)) 
    if len(A)!=0:
      data = Counter(A).most_common(1)[0][0] 
      #找出列表A中出现次数最多的元，如果有多个，则选择节点标号小的
      S=S|{data}   #将选出的节点加入种子集S 
      B=[B[i] for i in range(len(B)) if (data not in B[i])] 
      #删除选出的节点data所在的那些可达集
  return list(S)
def Sampling(G,k,ep,l): #生成RR集
  ep1=ep*2**0.5 #elpsilon'
  a=(l*math.log(n)+math.log(2))**0.5  #alpha(5)
  b=((1-math.e**(-1))*(math.log(comb(n,k))+l*math.log(n)+math.log(2)))**0.5  #beta(5)
  ep2=ep*a/((1-math.e**(-1))*a+b) #elpsilon1
  la1=(2+2*ep1/3)*(math.log(comb(n, k))+l*math.log(n)+math.log(math.log(n,2)))*n*ep1**(-2)
  #系数lambda' (9)
  la2=2*n*((1-math.e**(-1))*a+b)**2*ep**(-2) #系数lambda*(6)
  #print "a=",a,"b=",b,"la1=",la1,"la2=",la2  
  R=[]
  LB=1  
  ep1=ep*2**0.5
  for i in range(1,int(math.log(n1,2))):
    x=n1/2**i  
    th=la1/x  #theta_i
    while len(R)<=th:        
      t=random.randint(0,n1-1)  #随机选一个实体t
      T=set()
      for j in range(m1):   #依次选取实体t在每层对应的节点的RR的并集          
       S=[t+n1*j]  #t+n1*j为实体t在第j层对应的节点
       T1=runIC (Q, S)   #生成第j层对应点的RR集
       T=T|set(T1)     #每层对应点的RR集求并集为RRE集
      T=list(T)     #将实体t的RRE集转换格式为列表  
      R.append(T)    #把th个实体分别对应的th个RRE集放入列表R
    S=NodeSelection(R,k)
    if n*F(R,S)>=(1+ep1)*x:
      LB=n*F(R,S)/(1+ep1)
      break
  th=la2/LB
  while len(R)<th:
      t=random.randint(0,n1-1)  #随机选一个实体t
      T=set()
      for j in range(m1):   #依次选取实体t在每层对应的节点的RR的并集          
       S=[t+n1*j]  #t+n1*j为实体t在第j层对应的节点
       T1=runIC (Q, S)   #生成第j层对应点的RR集
       T=T|set(T1)     #每层对应点的RR集求并集为RRE集
      T=list(T)     #将实体t的RRE集转换格式为列表  
      R.append(T)    #把th个实体分别对应的th个RRE集放入列表R  
  return R    
            
def Inf(G,S): #蒙特卡罗方法算影响力
  num=0
  iternum=1000
  for i in range(iternum):
    T=runIC(G,S)
    T=len(set([T[i]%n1 for i in range(len(T))]))
    num+=T  
  num=num/iternum  
  return num



#公平性种子颜色分配问题
def Dic(t,n1):      #模n1同余的点（key）的概率值（value）叠加到一点（实体）t为字典
  s1=dict()
  for i in t:
    t1={key:value for key,value in t.items() if key%n1==i%n1}
    s= sum(v for v in t1.values())    
    s1[i%n1]=s
  return s1   #输出为字典，只保留点0到n1和他们的value
def runMIC (G, S1,S2):
    from copy import deepcopy
    from random import random
    T1,U1=set(),set()
    for i in S1:
      for j in range(m1):  
        T1=T1|{i%n1+j*n1}  
    for i in S2:
      for j in range(m1):  
        U1=U1|{i%n1+j*n1}
    T1=list(T1)   
    U1=list(U1)      #生成种子主体
    #print "T1,U1=",T1,U1
    T = deepcopy(S1) # copy already selected nodes
    U = deepcopy(S2)
    dict1 = {}
    i = 0
    while i < len(T):
        for v in G[T[i]]: # for neighbors of a selected node T[i]的邻点v
            if v not in T and v not in U1 and v not in T1: # if it wasn't selected yet，且非种子
                w = G[T[i]][v]['weight'] # count the number of edges between two nodes
                if random() <= w: # if at least one of edges propagate influence
                    #print T[i], 'influences', v                    
                    T.append(v)            #T是list文件        
                    d1 = {v: w} 
                    dict1.update(d1) 
                    #print "G[T[i]]:",G[T[i]] ,"w:",w,",1 - (1-p)**w:",1 - (1-p)**w
        i += 1
    dict2 = {}    
    j = 0
    while j < len(U):
        for v in G[U[j]]: # for neighbors of a selected node T[j]的邻点v
            if v not in U and v not in U1 and v not in T1 : # if it wasn't selected yet，且非种子
                w = G[U[j]][v]['weight'] # count the number of edges between two nodes
                if random() <= w: # 注意：不能在import random下使用，只在from random import random下
                    #print U[j], 'influences', v                    
                    U.append(v)                    
                    d2 = {v: w} 
                    dict2.update(d2) 
                    #print "G[T[i]]:",G[T[i]] ,"w:",w,",1 - (1-p)**w:",1 - (1-p)**w
        j += 1 
    #print "T,U=",T,U  
    T=list(set([c%n1 for c in T]))
    U=list(set([c%n1 for c in U]))
    #print "T0,U0=",T,U  
    #print "dict1,dict2=",dict1,dict2  
    dict1=Dic(dict1,n1)   
    dict2=Dic(dict2,n1)
    #print "dict11,dict22=",dict1,dict2  
    d1={key:value for key,value in dict1.items() if key in dict2}
    d2={key:value for key,value in dict2.items() if key in dict1}
    #print "d1,d2=",d1,d2
    for key in d1:
      t=random()  
      #print "random=",t
      #print "d1[key]/(d1[key]+d2[key])=",d1[key]/(d1[key]+d2[key])
      if t<= d1[key]/(d1[key]+d2[key]):         
        U.remove(key)
      else:
        T.remove(key)
    #print "T2,U2=",T,U 
    
    S1=[S1[l]%n1 for l in range(len(S1))]
    S2=[S2[l]%n1 for l in range(len(S2))]
    #print "S1,S2=",S1,S2
    S11,S22=set(S1),set(S2)
    #print "set S1,S2=",S11,S22
    for j in S11&S22:
      t=float(S1.count(j)) #整型转浮点型 
      u=float(S2.count(j)) 
      o=random()
      if o<=t/(t+u):
        U.remove(j)  
      else:
        T.remove(j) 
    return len(T),len(U)
     
def avgSize(G,S1,S2,iterations): #双重影响力结果
    avg = 0
    avg1 = 0
    for i in range(iterations):
        result=runMIC(G,S1,S2)   #元组数据            
        avg += float(result[0])/iterations           
        avg1 += float(result[1])/iterations    
    return avg,avg1,avg+avg1

def fair(d1,y,S):
  t9=time.time()
  ly=len(y)  #颜色的个数
  ls=len(S) #种子点个数
  d3 = d1.copy()       #颜色分配方法二
  lists=[[] for _ in range(ly)]
  f0=[0 for _ in range(ly)]  #初始化各个颜色的影响力
  for i in range(ls):
    a=min(d3, key=d1.get)   #d1中value最小的项
    f1=[f0[i]/y[i] for i in range(ly)] #计算各个颜色影响力和预算之比值（产出投入比）
    b=min(f1)         #f1中的最小值
    c=f1.index(b)     #f1最小值的index
    lists[c].append(a)  #将d1中value最小的项a添加到lists的第c项
    f0[c]+=d1[a]   #更新f0的第c项的值  
    del d3[a]
  f1=[f0[i]/y[i] for i in range(ly)] #f1最后一次更新
  print 'k=',ls,'t=',ly,'y=',y 
  print 'method2, the seed allocation is:',lists 
  #print 'respectively, the influences of colors are:',f1
  print'fair of method2=',min(f1)/max(f1) 
  t10=time.time()
  print "time of method2:", t10-t9  

  
  import random
  from random import randint
  f=0
  num=1000
  for j in range(num):
    d2 = d1.copy()
    lists=[[] for _ in range(ly)]
    f0=[0 for _ in range(ly)]  #初始化各个颜色的影响力
    for i in range(ls):      
      a=random.choice(d2.keys())   #随机取d1中一种子 
      c=randint(0, ly-1)   #随机取一种颜色
      lists[c].append(a)  #将种子a添加到颜色c
      f0[c]+=d2[a]   #更新颜色c的影响力
      f11=[f0[i]/y[i] for i in range(ly)] #计算各个颜色影响力和预算之比值（产出投入比）
      del d2[a]
    f+=min(f11)/max(f11)
  f=f/num  
  #print 'y=',y
  #print 'the seed allocation is:',lists
  #print 'respectively, the influences of colors are:',f1
  print'fair of method1(random)=', f 
  

  t13=time.time()
  d4 = d1.copy()     #方法3颜色分配
  f2=[1/float(y[i])  for i in range(ly)] #预算的倒数
  lists=[[] for _ in range(ly)]
  f0=[0 for _ in range(ly)]  #初始化各个颜色的影响力
  for i in range(ls):
    a=max(d4, key=d1.get)   #d1中value最大的项  
    b=min(f2)         #f2中的最小值
    c=f2.index(b)     #f2最小值的index,即颜色
    lists[c].append(a)  #将d1中value最大的项a添加到lists的第c项
    f0[c]+=d1[a]   #更新f0的第c项的值 
    del d4[a]
    f2[c]=f0[c]/y[c]
  print 'method3, the seed allocation is:',lists 
  print 'fair of method3=',min(f2)/max(f2) 
  #print 'respectively, the influences of colors are:',f2
  t14=time.time()
  print "time of method3:",t14-t13 
  
 #主程序（找种子节点）
for k in range(1,5):  
  print 'k=', k 
  t1=time.time()
  l=l*(1+math.log(2)/math.log(n))
  R=Sampling(G,k,ep,l)
  S=NodeSelection(R,k)
  print 'RREseeds are:',S  
  t2=time.time()
  print "RRE time:",t2-t1 
  infl=F(R,S)*n1
  print "the influence of RRE method:", infl 
  print "Monte Carlo influence of RRE method:", Inf(G,S) 
  #基于最大出度的方法
  t3=time.time()
  Seeds=[]
  G1=copy.deepcopy(G)
  for i in range(k):  
    Ds=dict(G1.out_degree())  #节点和对应出度的字典文件
    Mx=max(Ds.values())    #最大出度值
    t=max(Ds.iterkeys(), key=(lambda k: Ds[k]))  
    #最大出度的点,这里的k只是变量，和前文的k无关
    Seeds.append(t) #加入种子集
    G1.remove_node(t) #图G
  print 'Degseeds are:',Seeds 
  t4=time.time()
  print "time:",t4-t3 
  print 'influence of degree method:',Inf(G,Seeds) 
  t5=time.time()
  Ranseeds=random.sample(range(0, n-1), k)
  print 'ranseeds are:',Ranseeds 
  t6=time.time()
  print "time:",t6-t5 
  print 'influence of random method:',Inf(G,Ranseeds) 

  t7=time.time()
  d1=dict()
  print 'Seeds are:', S 
  ls=len(S) #种子点个数
  for i in S:    
    S1=[i]
    S2=[S[j] for j in range(ls) if S[j]!=i ]
    #print 'S1,S2,S=',S1,S2,S
    result=avgSize(G,S1,S2,1000)
    d1.update({i:result[0]})
    #print 'result=',result
    #print 'result[0]=',result[0]
    #print 'd1=',d1
  t8=time.time()
  print "time of computing single influence:",t8-t7    
  
  if k>1:  
    y=[1,1]  #多影响（颜色）的预算比例 作为输入2
    fair(d1,y,S)
    y=[1,2]  #多影响（颜色）的预算比例 作为输入2
    fair(d1,y,S)
    if k>2:
      y=[1,1,1]  #多影响（颜色）的预算比例 作为输入2
      fair(d1,y,S)
      y=[1,2,3]  #多影响（颜色）的预算比例 作为输入2
      fair(d1,y,S)
      if k>3:
        y=[1,1,1,1]  #多影响（颜色）的预算比例 作为输入2
        fair(d1,y,S)
        y=[1,2,3,4]  #多影响（颜色）的预算比例 作为输入2
        fair(d1,y,S)  
print 'done'   