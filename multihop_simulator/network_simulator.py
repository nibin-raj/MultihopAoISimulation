import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

from scipy.io import loadmat
from scipy import linalg as la

import logging
import itertools
from scipy.special import comb

from .constants import *
from .age_difference_scheduler import *
from .Index_scheduler import *


class Packet:
    def __init__(self, sourceid, destinationid, seqnum, generated_time):
        self.sourceid = sourceid
        self.destinationid = destinationid
        self.seqnum = seqnum
        self.generated_time = generated_time
        logging.info("Packet:from %u: to %u: seq %u: slot %u" % (sourceid, destinationid, seqnum, generated_time))
        
    def get_packetage(self, current_time):
        return current_time - self.generated_time + 1 # add +1; pkt generated in slot t can be served in t itself, but takes one slot time
    
class Node:
    def __init__(self, nodeid, totalnum_nodes, issource, destinationid, scheduler_id):
        self.buffer = []
        for i in range(totalnum_nodes):
            self.buffer.append([])
    
        self.nodeid = nodeid
        self.scheduler_id=scheduler_id
        self.packet_seqnum = 0
        self.issource = issource
        self.destinationid = destinationid
        logging.info("Node:Created:ID %u: IsSrc %u: DestId %u" % (nodeid, issource, destinationid))

        self.last_removed_packet = [NOPACKET] * totalnum_nodes

                
    def generate_packet(self, current_time):
        if self.issource:
            packet = Packet(self.nodeid, self.destinationid, self.packet_seqnum, current_time)
            self.buffer[self.nodeid].append(packet)
            self.packet_seqnum += 1
            logging.info("Node:Packet g+:ID %u: PacketSeq %u: Bufferlength %u" % (self.nodeid, packet.seqnum, len(self.buffer[self.nodeid])))
            
    def add_packet(self, packet):
        # assuming tree like topologies
        # self.buffer[packet.sourceid].append(packet)
        self.buffer[packet.sourceid] = [packet]
        logging.info("Node:Packet +:ID %u: PacketSeq %u: Bufferlength %u" % (self.nodeid, packet.seqnum, len(self.buffer[packet.sourceid])))
        
    def add_packet_to_root(self, packet):
        self.buffer[packet.sourceid] = [packet]
        logging.info("Node:Packet +:ID %u: PacketSeq %u: Bufferlength %u" % (self.nodeid, packet.seqnum, len(self.buffer[packet.sourceid])))
        
        
    def remove_packet_from_hol(self, sourceid):
        if len(self.buffer[sourceid]) > 0:
            packet = self.buffer[sourceid].pop(0)
            logging.info("Node:Packet -:ID %u: PacketSeq %u: Bufferlength %u" % (self.nodeid, packet.seqnum, len(self.buffer[packet.sourceid])))
            self.last_removed_packet[sourceid] = packet
            return packet
        else:
            return NOPACKET
        
    def get_latest_received_packet(self, sourceid):
        if len(self.buffer[sourceid]) > 0:
            return self.buffer[sourceid][-1]
        else:
            return NOPACKET
        
    def logmeasurements_oneslot(self, t):
        ages = []
        for si, b in enumerate(self.buffer):
            pkt = self.get_latest_received_packet(si)
            if not pkt == NOPACKET:
                age = pkt.get_packetage(t)
            else:
                if not self.last_removed_packet[si] == NOPACKET:
                    age = self.last_removed_packet[si].get_packetage(t)
                else:
                    age = 0
            ages.append(age)
        ages = ",".join([str(i) for i in ages])
        logging.info("AgeMeasurement:%u,%u,%u,%s" % (t, self.scheduler_id, self.nodeid, ages))
        
    

class Network:
    def __init__(self, totalnum_nodes, link_list, source_list, commissioned_nodes, network_type, scheduler_id, interference_model, interference_k_hop):
        self.totalnum_nodes = totalnum_nodes
        self.link_list = link_list
        self.source_list = source_list
        self.commissioned_nodes = commissioned_nodes
        
        self.network_type = network_type
        self.interference_model = interference_model
        self.interference_k_hop = interference_k_hop
        
        
        self.G, self.G_up, self.G_down, self.leaf_nodes, self.line_graphs, self.subtree_roots = self.make_graph()
        self.A = self.get_activation_vectors()
        
        for n in self.G_up.nodes:
            issource = False
            if n in self.source_list:
                issource = True
                
            node = Node(n, self.totalnum_nodes, issource, ROOTNODE_ID, scheduler_id)
            self.G_up.nodes[n]["Node"] = node
        
    def make_graph(self):
        G_up = nx.DiGraph()
        G_down = nx.DiGraph()
        G = nx.Graph()
        root = 0
        for li, l in enumerate(self.link_list):
            G_up.add_edge(l[1], l[0])
            G_down.add_edge(l[0], l[1])
            G.add_edge(l[1],l[0])
            
        leaf_nodes = []
        for n, d in G_up.in_degree:
            if d == 0:
                leaf_nodes.append(n)
        
        line_graphs = []
        for n in sorted(leaf_nodes):
            line_graphs.append(nx.shortest_path(G_up, n, 0))
        
        subtree_roots = [n for n in nx.neighbors(G_down, 0)]

        return G, G_up, G_down, leaf_nodes, line_graphs, subtree_roots
    
    def logmeasurements_oneslot(self, t):
        for n in self.G_up.nodes:
            self.G_up.nodes[n]["Node"].logmeasurements_oneslot(t)

    def update_node_state(self, state):
        state_dict = {"R":0, "T":1, "I":2}
        next_state = ["T","I","R"]
        return next_state[state_dict[state]]
        
    def get_activation_vectors_only_edge(self):
        radio_state = ["R","T","I"]
        Mat = []
        for b in sorted(self.line_graphs):
            for n in sorted(b):
                if n == 0:
                    continue
                self.G_up.nodes[n]["State"] = radio_state[nx.shortest_path_length(self.G_up, n, 0) % 3]

        for aa in range(3):
            for b in sorted(self.line_graphs):
                for n in sorted(b):
                    if n == 0:
                        continue
                    if self.G_up.nodes[n]["State"] == "T":
                        Mat.append(1)
                    else:
                        Mat.append(0)
                    self.G_up.nodes[n]["State"] = self.update_node_state(self.G_up.nodes[n]["State"])

        for n in self.G_up.nodes():
            if n == 0:
                continue
            del self.G_up.nodes[n]["State"]

        Mat = np.array(Mat)
        M_mat0 = Mat.reshape(3, len(self.G_up.nodes) - 1)
        return M_mat0
    
    def get_activation_vectors_singleline(self):
        radio_state = ["R","T","I"]
        Mat = []
        for b in sorted(self.line_graphs):
            for n in sorted(b):
                if n == 0:
                    continue
                self.G_up.nodes[n]["State"] = radio_state[nx.shortest_path_length(self.G_up, n, 0) % 3]

        for aa in range(3):
            for b in sorted(self.line_graphs):
                for n in sorted(b):
                    if n == 0:
                        continue
                    if self.G_up.nodes[n]["State"] == "T":
                        Mat.append(1)
                    else:
                        Mat.append(0)
                    self.G_up.nodes[n]["State"] = self.update_node_state(self.G_up.nodes[n]["State"])

        for n in self.G_up.nodes():
            if n == 0:
                continue
            del self.G_up.nodes[n]["State"]

        Mat = np.array(Mat)
        M_mat0 = Mat.reshape(3, len(self.G_up.nodes) - 1)

        awithsrc = []
        
        for i in range(M_mat0.shape[0]):
            a = M_mat0[i]
            l = len(np.where(a)[0])
            if l == 0:
                continue
            for j in itertools.product(self.source_list, repeat = l):
                ta = a.copy()
                ta[np.where(a)[0]] = j
                awithsrc.append(list(ta))
        return awithsrc
    
    def get_edgenode_graph(self):
        H = nx.Graph()
        for u,v in self.G_up.edges:
            H.add_node((u,v))
            
        hnodes = list(H.nodes)
        l_hnodes = len(hnodes)
        for i in range(l_hnodes):
            for j in range(i + 1, l_hnodes):
                s1 = set(hnodes[i])
                s2 = set(hnodes[j])
                if len(s1.intersection(s2)) > 0:
                    H.add_edge(hnodes[i], hnodes[j])
        return H
    
    def get_edgenode_graph_ind(self):
        H = nx.Graph()
        for u,v in self.G_up.edges:
            H.add_node((u,v))
            
        hnodes = list(H.nodes)
        l_hnodes = len(hnodes)
        for i in range(l_hnodes):
            for j in range(i + 1, l_hnodes):
                s1 = set(hnodes[i])
                s2 = set(hnodes[j])
                if 0 in s1 and 0 in s2:
                    continue
                if len(s1.intersection(s2)) > 0:
                    H.add_edge(hnodes[i], hnodes[j])
        return H
    
        
    def get_interference_matrix(self):
        num_links = len(self.link_list)
        links = list(self.G_up.edges)
        
        if self.interference_model == COMPLETE_INT:
            interference_matrix = np.ones((num_links, num_links))
            interference_matrix = interference_matrix - np.eye(num_links)
            return interference_matrix
        
        if self.interference_model == K_HOP_INT:
            H = self.get_edgenode_graph()
            interference_matrix = np.zeros((num_links, num_links))
            link_to_index_map = {}
            for li, l in enumerate(links):
                link_to_index_map[l] = li

            for li, l in enumerate(links):
                khop_graph = nx.ego_graph(H, l, undirected = True, radius = self.interference_k_hop)
                khop_neighbors = list(khop_graph.nodes)
                for n in khop_neighbors:
                    interference_matrix[li, link_to_index_map[n]] = 1
            interference_matrix = interference_matrix - np.eye(num_links)
            return interference_matrix
        
        if self.interference_model == K_HOP_INT_IND:
#             print('K_HOP_INT_IND')
            H = self.get_edgenode_graph_ind()
            interference_matrix = np.zeros((num_links, num_links))
            link_to_index_map = {}
            for li, l in enumerate(links):
                link_to_index_map[l] = li

            for li, l in enumerate(links):
                khop_graph = nx.ego_graph(H, l, undirected = True, radius = self.interference_k_hop)
                khop_neighbors = list(khop_graph.nodes)
                for n in khop_neighbors:
                    interference_matrix[li, link_to_index_map[n]] = 1
            interference_matrix = interference_matrix - np.eye(num_links)
            for b in self.subtree_roots:
                l = (b,0)
                for b_dash in self.subtree_roots:
                    l_dash = (b_dash,0)
                    interference_matrix[ link_to_index_map[l], link_to_index_map[l_dash]] = 1
            return interference_matrix
                
        if self.interference_model == GENERAL_INT:
            interference_matrix = np.zeros((num_links, num_links))
            return interference_matrix

            
    def is_independent_set(self, graph, nodes):
        for node1 in nodes:
            for node2 in nodes:
                if node1 != node2 and graph.has_edge(node1, node2):
                    return False
        return True

    def find_all_independent_sets(self, graph):
        independent_sets = []
        for nodes in graph.nodes():
            if self.is_independent_set(graph, [nodes]):
                independent_sets.append([nodes])

        for size in range(2, len(graph.nodes()) + 1):
            for nodes in combinations(graph.nodes(), size):
                if self.is_independent_set(graph, nodes):
                    independent_sets.append(list(nodes))
        return independent_sets


    def get_activation_vectors(self):
        interference_matrix = self.get_interference_matrix()
        # print("IMatrix", interference_matrix)
        r, c = np.where(interference_matrix)
        interference_graph = nx.Graph(zip(r,c))
        
        if self.network_type == multisource:
#             print('Multisource')
            interference_graph_complement = nx.complement(interference_graph)
            cliques = nx.find_cliques(interference_graph_complement)
    #         print("Cliques", cliques)
            links = list(self.G_up.edges)
            activation_vector_list = []

            for c in cliques:
                activation_vector = []
                for li in c:
                    link = links[li]
                    activation_vector.append(link)
                activation_vector_list.append(activation_vector)

            
            activation_vector_list_with_src = []
            for av in activation_vector_list:
                lav = len(av)
#                 print(av, lav)
                for j in itertools.product(self.source_list, repeat = lav):
                    avs = []
                    for l in range(lav):
                        if  j[l] in list(nx.ancestors(self.G_up, av[l][0])) or (j[l] == av[l][0]):
                            avs.append([j[l], av[l][1], av[l][0]]) 
                        else:
                            continue
                        if len(avs)==lav:
                            activation_vector_list_with_src.append(avs)
                            
                    
        if self.network_type == singlesource:
#             print('Singlesource')
            independent_sets = self.find_all_independent_sets(interference_graph)
            links = list(self.G_up.edges)
            
            activation_vector_list = []
            for c in independent_sets:
                activation_vector = []
                for li in c:
                    link = links[li]
                    activation_vector.append(link)
                activation_vector_list.append(activation_vector)

            activation_vector_list_with_src = []
            for av in activation_vector_list:
                lav = len(av)
                for j in itertools.product(self.source_list, repeat = lav):
                    avs = []
                    for l in range(lav):
                        avs.append([j[l], av[l][1], av[l][0]])
                    activation_vector_list_with_src.append(avs)
                
        return activation_vector_list_with_src
    
        
    def get_link_s_ht(self, activation_vector):
        ht = []
        for li, l in enumerate(activation_vector):
            if l > 0:
                ht.append([l, li, li + 1])
        return ht
    
    def isfeasible(self, activation_vector):
        # sht = self.get_link_s_ht(activation_vector)
        sht = activation_vector
        for (s, h, t) in sht:
            if (s != t) and len(self.G_up.nodes[t]["Node"].buffer[s]) == 0:
                return False
        return True
    
    def pack_age(self, source, current_time):
        last_received_packet = self.G_up.nodes[ROOTNODE_ID]["Node"].get_latest_received_packet(source)
        if not last_received_packet == NOPACKET:
            packetageis = last_received_packet.get_packetage(current_time)
            return packetageis
        else:
            return current_time + 1

        

class NetworkSimulator:
    def __init__(self, totalnode_num, link_list, source_list, commissioned_nodes, psv, network_type, scheduler_id, scheduler, interference_model, interference_k_hop):
        self.totalnode_num = totalnode_num
        self.link_list = link_list
        self.source_list = source_list
        self.commissioned_nodes = commissioned_nodes
        self.ps_values = psv
        self.network_type = network_type
        
        self.scheduler = scheduler
        self.scheduler_id = scheduler_id
        self.network = Network(self.totalnode_num, self.link_list, self.source_list, self.commissioned_nodes, self.network_type,self.scheduler_id, interference_model, interference_k_hop)
        self.initialize_scheduler()
        
        self.packetagedata = []
        for i in range(self.network.totalnum_nodes):
            self.packetagedata.append([])
        
    def initialize_scheduler(self):
        self.scheduler.network = self.network
        if self.scheduler_id == SCHEDULER_RANDOM:
            t = len(self.network.A)
            d = [1.0/t] * t
            self.scheduler.activation_vector_distribution = d
            self.scheduler.setup_nodedata()
        if self.scheduler_id == SCHEDULER_AGEDEBT:
            self.scheduler.setup_nodedata()
        if self.scheduler_id == SCHEDULER_AGEDIFF:
            self.scheduler.setup_nodedata()
        
        
    def simulate_oneslot(self, t):
        logging.info("Simulator:Slot %u" % (t))
        source_generatedpacket = []
        activation_vector_index = NOT_SELECTED
        if self.scheduler_id == SCHEDULER_RANDOM:
            activation_vector_index, source_generatedpacket = self.scheduler.get_activation_vector_slot()
        if self.scheduler_id == SCHEDULER_AGEDEBT:
            activation_vector_index, source_generatedpacket = self.scheduler.get_activation_vector_slot()
        if self.scheduler_id == SCHEDULER_AGEDIFF:
            activation_vector_index, source_generatedpacket, randomvalue = self.scheduler.get_activation_vector_slot(self.ps_values)
        
        logging.info("Simulator:Slot %u: Source Gen Packets %s" % (t, source_generatedpacket))
        
        for s in source_generatedpacket:
            self.network.G_up.nodes[s]["Node"].generate_packet(t) 
            
        if not activation_vector_index == NOT_SELECTED:
            activation_vector = self.network.A[activation_vector_index]
            logging.info("Simulator:Slot %u: Act. vector %s" % (t, activation_vector))
            delivery_t = t 
            # shts = self.network.get_link_s_ht(activation_vector)
            shts = activation_vector
            for (s, h, t) in shts:
                packet = self.network.G_up.nodes[t]["Node"].remove_packet_from_hol(s)
                if not packet == NOPACKET:
                    if not h == ROOTNODE_ID:
                        if self.scheduler_id ==SCHEDULER_AGEDIFF:
                            if randomvalue < self.ps_values[(h,t)]:
                                self.network.G_up.nodes[h]["Node"].add_packet(packet)
                        else:
                            self.network.G_up.nodes[h]["Node"].add_packet(packet)
                            
                    else:
                        if self.scheduler_id ==SCHEDULER_AGEDIFF:
                            if randomvalue < self.ps_values[(h,t)]:
                                self.network.G_up.nodes[h]["Node"].add_packet_to_root(packet)
                        else:
                            self.network.G_up.nodes[h]["Node"].add_packet_to_root(packet)

                    
        for s in self.network.source_list:
            find_packetage = self.network.pack_age(s, delivery_t)
            self.packetagedata[s].append(find_packetage)
                    
    def logmeasurements_oneslot(self, t):
        self.network.logmeasurements_oneslot(t)

    def simulate(self, max_slots):
        for t in range(max_slots):
            self.simulate_oneslot(t)
            self.logmeasurements_oneslot(t)
