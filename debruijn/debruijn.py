#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

# ~/.local/bin/

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
import textwrap
matplotlib.use("Agg")

__author__ = "DELHOMME Jean"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["DELHOMME Jean"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "DELHOMME Jean"
__email__ = "jean.delhomme@live.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):

    with open (fastq_file, "r") as file:
        for ligne in file:
            yield next(file).strip()
            next(file)
            next(file)
            


def cut_kmer(read, kmer_size):

    for i in range(len(read) - kmer_size + 1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):

    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
        
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
        
    return (kmer_dict)



def build_graph(kmer_dict):

    digraph = nx.DiGraph()
    for kmer, count in kmer_dict.items():
        digraph.add_edge(kmer[:-1], kmer[1:], weight=count)
    
    return digraph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):

    for node in path_list:

        if delete_entry_node == True and delete_sink_node == True:
            graph.remove_nodes_from(node)
            #tous les noeuds d un chemin sont supprimes

        elif delete_entry_node == True:
            graph.remove_nodes_from(node[:-1])
            #tous les noeuds d un chemin sont supprime sauf le dernier

        elif delete_sink_node == True:
            graph.remove_nodes_from(node[1:])
            #tous les noeuds d un chemin sont supprmies sauf le premier

        elif delete_entry_node == False and delete_sink_node == False:
            graph.remove_nodes_from(node[1:-1])
            #tous les noeuds d un chemin sont supprimes saud le premier et le dernier

    return graph

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):

    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    
    starting_node_list = []
    
    for node in graph.nodes():
         if not list(graph.predecessors(node)):
            starting_node_list.append(node)

    return(starting_node_list)

def get_sink_nodes(graph):

    sink_node_list = []
    
    for node in graph.nodes():
         if not list(graph.successors(node)):
            sink_node_list.append(node)

    return(sink_node_list)

def get_contigs(graph, starting_nodes, ending_nodes):

    contig = []

    for node_start in starting_nodes:
        for node_end in ending_nodes:

            if nx.has_path(graph, node_start, node_end) == True:
                for path in nx.all_simple_paths(graph, node_start, node_end):
                    seq = path[0]
                    for node in path[1:]:
                        seq += node[-1]
                    contig.append(tuple((seq, len(seq))))

    return contig

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as file:
        for i, (contig, length) in enumerate(contigs_list):
            print(contig)
            #file.write(f">contig_{i} len={length}\n{fill(contig, width=80)}\n")
            file.write(">contig_{} len={}\n{}\n".format(i, length, textwrap.fill(contig, width = 80)))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    read_generator = read_fastq(args.fastq_file)

    for read in read_generator:
        cut_kmer(read, args.kmer_size)
    build_kmer_dict(args.fastq_file, args.kmer_size)

    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))

    nx.draw(graph)
    plt.show()

    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)

    save_contigs(contigs_list, args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
    
    

