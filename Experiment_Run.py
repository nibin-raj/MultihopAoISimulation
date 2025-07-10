from multihop_simulator import *

link_list = [(0, 1), (0, 9), (1, 4), (1, 8), (5, 6), (7, 5), (9, 7), (9, 10), (10, 2), (10, 3)]

# Network_Type = singlesource 
Network_Type = multisource 
interference_k_hop = 2
iterno = 1000

G_up = nx.DiGraph()
G_down = nx.DiGraph()
root = 0
for l in link_list:
    G_up.add_edge(l[1], l[0])
    G_down.add_edge(l[0], l[1])
            
        
commissioned_nodes = {}
for node in G_down.nodes():
    if node ==0: continue
    ancestors = nx.ancestors(G_down, node)
    ancestors_list = sorted(list(ancestors))
    ancestors_list = [x for x in ancestors_list if x != 0]
    commissioned_nodes[node] = ancestors_list
            
source_list = []
for node in G_up.nodes():
    if node ==0: continue 
    source_list.append(node)
    
totalnum_nodes = len(link_list) + 1

ps_values = {}
for sou in link_list:
    head_,tail_ = sou[0],sou[1]
    ps_values[(head_,tail_)] = 0.6  #Homogeneous Links
source_ps = list(ps_values.values())

source_hs = [nx.shortest_path_length(G_up, source=node, target=0) for node in source_list]

#### Index Polciy
num_sources = len(source_list)
sources_index = [Source(source_hs[i],source_ps[i]) for i in range(num_sources)]
for t in range(iterno): 
    indices = np.zeros((num_sources, 1)) 
    for s in range(num_sources):  
        indices[s] = sources_index[s].compute_index(source_ps[s])
    scheduled_source = np.argmax(indices)

    for s in range(num_sources):
        if s == scheduled_source: 
            sources_index[s].step(True, source_ps[s])    
        else:
            sources_index[s].step(False, source_ps[s])
scheduler_aoi = [np.mean(source.age_track) for source in sources_index] 
AoI_Indexscheduler = np.mean(scheduler_aoi)
print(AoI_Indexscheduler)

#### Age-Difference
agediffscheduler = AgeDifferenceScheduler()
simulator = NetworkSimulator(totalnum_nodes, link_list, source_list, commissioned_nodes, ps_values, Network_Type, SCHEDULER_AGEDIFF, agediffscheduler, COMPLETE_INT, interference_k_hop)
simulator.simulate(iterno)
av1 =[]
for n in source_list:
    av1.append(np.sum(agediffscheduler.agedata[n])/len(agediffscheduler.agedata[n]))
Av_agediff =  np.mean(av1)
print(Av_agediff)    
    
    
    