import os
import csv
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas
import pandas as pd
import logging
from .startup import HGNC_data_path, UniProt_data_path
from typing import Dict, Tuple, Iterable, Optional

logger = logging.getLogger('network')



class Network():
    def __init__(self,
                 nodes: dict,
                 graph: nx.Graph,
                 ppi_file: Optional[str],
                 node_path: Optional[str],
                 edge_path: Optional[str]):
        self.nodes = nodes
        self.graph = graph
        self.ppi_file = ppi_file
        self.node_path = node_path
        self.edge_path = edge_path

        # initialize the methods, user can put either PPIs file or Node/Edge list, and update the nodes dictionary.
        if self.ppi_file and not self.node_path and not self.edge_path:
            self.generate_node_dict(self.read_ppis(self.ppi_file))

        elif not self.ppi_file and self.node_path and self.edge_path:
            self.original_dict =self.node_label(self.relations(self.node_path))

        # Raise ERROR when three files are there.
        elif self.ppi_file and self.node_path and self.edge_path:
            raise IOError('Cannot import three files!')

    def read_ppis(self, ppi_file: str) -> Iterable[Tuple[str]]:
        """ Read original PPI file
        Parameters
        ----------
        ppi_file: str
                 .csv PPI file
        Returns
        -------
        Iterable[Tuple[str]]
                The relationship of nodes and type of interaction.
        """
        f = open(ppi_file)
        content = f.readlines()
        self.rels = [tuple(x.strip().split(',')) for x in content[1:]]
        return self.rels

    def generate_node_dict(self, relations: Iterable[Tuple[str]]) -> Dict[str, str]:
        """
           From the relation tuple to node dictionary, where key is the identifier and value is the HGNC symbol
        Parameters
        ----------
        relations: Iterable[Tuple[str]]
                  from read_ppis function.
        Returns
        -------
        dict: A dictionary contains identifier as key and HGNC symbol as value
        """
        nodes_ = set()
        for rel in relations:
            nodes_.add(rel[0])
            nodes_.add(rel[2])
        self.nodes = {index: symbol for index, symbol in enumerate(nodes_, 1)}
        return self.nodes

    # Generate node list
    def write_node_list(self, node_path: str) -> None:
        """
            Write the identifier and symbol to .tsv file.
        Parameters
        ----------
        node_path: str
                 .tsv file path of node list.
        Returns
        -------
        None
        """
        with open(node_path, 'wt') as outfile_1:
            tsv_writer_node = csv.writer(outfile_1, delimiter='\t')
            for identifier, symbol in self.nodes.items():
                tsv_writer_node.writerow([identifier, symbol])

    def write_edge_list(self, edge_path: str) -> None:
        """
            Write interaction information to edge list file
        Parameters
        ----------
        edge_path: str
                  .tsv file of edge list
        Returns
        -------
        None
        """
        self.nodes_new = {symbol: identifier for identifier, symbol in self.nodes.items()}
        with open(edge_path, 'wt') as outfile_2:
            tsv_writer_edge = csv.writer(outfile_2, delimiter = '\t')
            for rel in self.rels:
                out_, in_ = self.nodes_new[rel[0]], self.nodes_new[rel[2]]
                interaction = rel[1].replace(" ", "_")
                tsv_writer_edge.writerow([out_, in_, interaction])

    def relations(self, node_path: str) -> Iterable[Tuple[str]]:
        """
        If giving node_list or edge list, then generating the corresponding relations from these files.
        Parameters
        ----------
        node_path: str
                 .tsv file of node/edge list
        Returns
        -------
        Iterable[Tuple[str]]
        The relations of identifier and HGNC symbol
        """
        f = open(node_path)
        content = f.readlines()
        self.new_rels = [tuple(x.strip().split('\t')) for x in content]
        return self.new_rels

    def node_label(self, relations: Iterable[Tuple[str]]) -> dict:
        """
        From relationship Iterable, to the node dictionary, where identifier is key and HGNC symbol is value.
        Parameters
        ----------
        relations: Iterable[Tuple[str]]
                  The result from relations function.
        Returns
        -------
        dict
        A dictionary contains identifier and HGNC symbol.
        """
        for identifier, symbol in relations:
            self.nodes[int(identifier)] = symbol
        return self.nodes

    def import_graph(self, edge_path: str) -> nx.Graph:
        """
        To import nx.graph from edge file and use node dictionary to replace HGNC symbol with identifier.
        Parameters
        ----------
        edge_path: str
                  The .tsv edge file of PPIs.
        Returns
        -------
        nx.Graph
        """
        self.Data = pd.read_csv(edge_path, sep='\t', names=['node_1', 'node_2', 'metadata'])
        self.new_data = self.Data.replace(self.nodes)

        self.graph = nx.from_pandas_edgelist(self.new_data, source='node_1', target='node_2', edge_attr='metadata')
        return self.graph

    def check_output(self, graph_output: str) -> None:
        """
        To check if output extension is available.
        Parameters
        ----------
        graph_output: str
                     The graph file path.
        Returns
        -------
        None
        """
        accept_extension = ['pdf', 'png', 'jpg', 'svg']
        base_name = os.path.basename(graph_output)
        output_extension = base_name.split(".")[1]
        if output_extension not in accept_extension:
            raise ValueError("Graph extension is wrong, must use 'pdf', 'svg','png', 'jpg'.")

    def generate_graph_network(self, graph_output: str, print_edge_label: bool = False) -> None:
        """
        Use nx.Graph to generate graph PPIs network.
        Parameters
        ----------
        graph_output: str
                     The correct output file of network.
        print_edge_label: bool
                        If true, will print the edge label (type of protein interaction).
        Returns
        -------
        None
        """
        self.check_output(graph_output)
        node_colors, edge_colors = 'red', 'black'
        plt.figure(figsize=(18, 18))
        self.graph.pos = nx.spring_layout(self.graph, k = 0.06)
        nx.draw_networkx(self.graph, pos = self.graph.pos, with_labels=True,
                         node_color = node_colors,
                         edge_color = edge_colors, alpha = 0.3)

        # To check if need to print edge label
        if print_edge_label:
           edge_label = {(node1, node2): type['metadata'] for (node1, node2, type) in self.graph.edges(data=True)}
           nx.draw_networkx_edge_labels(self.graph, pos=self.graph.pos, edge_labels=edge_label)
        plt.savefig(graph_output)


class Analyzer(Network):

    def __init__(self,nodes: dict,
                 graph: nx.Graph,
                 ppi_file: Optional[str],
                 node_path: Optional[str],
                 edge_path: Optional[str]):
        super().__init__(nodes, graph, ppi_file, node_path, edge_path)
        self.short_path = []
        self.enrich_identifier_info = defaultdict(dict)

    def shortest_path(self, source: str, target: str, print_option: bool = False) -> Optional[list]:
        """
        Get shortest path between two given nodes.
        Parameters
        ----------
        source: str
               One Node
        target: str
               The other node
        print_option: bool
                     if ture, will print the shortest path in terminal.
        Returns
        -------
        Optional[list]
        The possible shortest path between two nodes.
        """
        self.paths = nx.all_shortest_paths(self.graph, source, target)
        self.short_path = [path for path in self.paths]

        if print_option:
            for elem in self.short_path:
                path_string = " -> ".join(elem)
                print('START : ***' + path_string + "  ***STOP")
        return self.short_path

    def color_path(self, path_nodes: list) -> tuple:
        """
        Generate different color for the path.
        Parameters
        ----------
        path_nodes: list
                   The list of nodes
        Returns
        -------
        tuple
        The customized color of nodes and edges
        """
        node_set = set()
        edge_set = set()
        for single_path in path_nodes:
            node_set.update(single_path)
            for i in range(len(single_path) - 1):
                edge_set.add((single_path[i], single_path[i + 1]))
        node_color_design = ['blue' if nodes in node_set else 'red' for nodes in self.graph.nodes()]
        edge_color_design = ['purple' if edges in edge_set else 'black' for edges in self.graph.edges()]
        return node_color_design, edge_color_design

    def generate_graph_network(self, graph_output: str, print_edge_label: bool = False) -> None:
        """
        Overwrite generate_graph_network function from the parent class,
        if there is short path, then color the path and nodes.
        Parameters
        ----------
        graph_output: str
                     The network output with shortest path.
        print_edge_label: bool
                     If ture, will print the edge label which illustrates the PPIs.
        Returns
        -------
        None
        """
        self.check_output(graph_output)  # check the extension of output file.
        # customize the color of edges for shortest paths.
        if self.short_path:
            node_colors, edge_colors = self.color_path(self.short_path)
        else:
            node_colors, edge_colors = 'red', 'black'

        plt.figure(figsize = (18, 18))
        self.graph.pos = nx.spring_layout(self.graph, k = 0.06)
        nx.draw_networkx(self.graph, pos=self.graph.pos,
                         with_labels=True,
                         node_color = node_colors,
                         edge_color = edge_colors)
        if print_edge_label:
            edge_label = {(node1, node2): type['metadata'] for (node1, node2, type) in self.graph.edges(data=True)}
            nx.draw_networkx_edge_labels(self.graph, pos=self.graph.pos, edge_labels=edge_label)
        plt.savefig(graph_output)

    def enrich_gather_identifier(self, node_path:str) -> dict:
        """
        This funciton is used to request HGNC symbol in the node_path file and extract info to a nested dictionary,
        where key is the HGNC symbol, and value is the dictionary of all HGNC id, emsembl id and uniprot id.
        -------
        parameter: node_path: str
                  original node_list path without metadata
        -------
        return: dict
               nested dictionary. key is the HGNC symbol, and values is a dictionary with its meta data.
        """
        from .Utils import Profiler
        from tqdm import tqdm
        f = open(node_path)
        content = f.readlines()
        rels = [tuple(x.strip().split('\t')) for x in content]
        symbol_list = [symbol for identifier, symbol in rels]

        # cache file : which is not present in HGNC and UniProt folder.
        # hgnc file
        for elem in tqdm(symbol_list):
            import time
            time.sleep(0.1)
            hgnc_path_root = HGNC_data_path
            hgnc_file_name = f"/{elem}.json"
            hgnc_path = hgnc_path_root + hgnc_file_name
            if not os.path.isfile(hgnc_path):
                Profiler(elem).request()
        # uniprot file
        for elem in tqdm(symbol_list):
            time.sleep(0.1)
            if os.path.isfile(HGNC_data_path + f"/{elem}.json"):
                self.enrich_identifier_info[f"hgnc_{elem}"] = Profiler(elem).extract()
                for acc_num in Profiler(elem).extract()["UniProt ID"]:
                    if acc_num is not None:
                        uniprot_path_root = UniProt_data_path
                        uniprot_file_name = f"/{acc_num}.fasta"
                        uniprot_path = uniprot_path_root + uniprot_file_name
                        if not os.path.isfile(uniprot_path):
                            Profiler(elem).request_uniprot(acc_num)
                        self.enrich_identifier_info[f"uniprot_{acc_num}"] = Profiler(elem).extract_uniprot(acc_num)
                    else:
                        self.enrich_identifier_info[f"uniprot_{acc_num}"] = {None}
        return

    def enrich_generate_node_dict(self, query_from_sql: bool = False) -> dict:
        """
        A function used to get node dictionary of enriched info( DNA, RNA and Protein)
        Returns
        -------
        dict
        A dictionary of enriched info.
        """
        self.nodes = defaultdict(dict)
        # A set contains unique node.
        nodes_ = set()
        for rel in self.rels:
            nodes_.add(rel[0])
            nodes_.add(rel[2])
        # A nested dictionary contains identifier as key, and in inner dict, the HGNC symbol as key.
        identifier = 1
        for symbol in nodes_:
            if query_from_sql:
                from .models import query_data
                if query_data(symbol):
                    gene_symbol = list(query_data(symbol).keys())[0]
                    value_list = list(query_data(symbol).values())[0]
                    ensembl = value_list[1]
                    acc_num = value_list[2]
                    tax_num = value_list[3]
                    self.nodes[identifier][symbol] = 'DNA'+" "+ f"HGNC:{gene_symbol} /"+ f"Ensembl:{ensembl}"
                    identifier += 1
                    self.nodes[identifier][symbol] = 'RNA'
                    identifier += 1
                    self.nodes[identifier][symbol] = 'Protein' + " " + f"{acc_num} /"+ f"Taxonomy:{tax_num}"
                    identifier += 1

            else:
            # if request gene and protein info from database:
                if self.enrich_identifier_info:
                    self.nodes[identifier][symbol] = 'DNA'+" "+ str(self.enrich_identifier_info[f"hgnc_{symbol}"]["HGNC ID"][0]) +"/"+ str(self.enrich_identifier_info[f"hgnc_{symbol}"]["Ensembl Gene ID"][0])
                    identifier += 1
                    self.nodes[identifier][symbol] = 'RNA'
                    identifier += 1
                    self.nodes[identifier][symbol] = 'Protein' + " " + str(self.enrich_identifier_info[f"hgnc_{symbol}"]["UniProt ID"][0])
                    identifier += 1
            # no gene and protein info:
                else:
                    self.nodes[identifier][symbol] = 'DNA'
                    identifier += 1
                    self.nodes[identifier][symbol] = 'RNA'
                    identifier += 1
                    self.nodes[identifier][symbol] = 'Protein'
                    identifier += 1


    def enrich_write_node_list(self, node_path: str) -> None:
        """
        Write node with their corresponding info( DNA, RNA, Protein) to the node_list file.
        Parameters
        ----------
        node_path: str
                  .tsv file with identifier and their info (DNA, RNA, Protein)
        Returns
        -------
        None
        """
        with open(node_path, 'wt') as outfile_1:
            tsv_writer_node = csv.writer(outfile_1, delimiter='\t')
            for identifier, info in self.nodes.items():
                for name, type in info.items():
                    tsv_writer_node.writerow([identifier, name, type])

    # from PPIs file:
    def enrich_write_edge_list(self, relation: Iterable[Tuple[str]], edge_path: str, query_from_sql: bool = False) -> None:
        """
        Write the relationship between two nodes to edge_list file.

        Parameters
        ----------
        relation: Iterable[Tuple[str]]
                 The relationship between two nodes
        edge_path: str
                 The output file of edge_list.
        Returns
        -------
        None
        """
        with open(edge_path, 'wt') as outfile_2:
            tsv_writer_edge = csv.writer(outfile_2, delimiter='\t')
            for n in range(len(relation) - 1):
                if relation[n][2][:3] == 'DNA':
                    tsv_writer_edge.writerow([relation[n][0], relation[n + 1][0], 'transcribed'])
                elif relation[n][2] == 'RNA':
                    tsv_writer_edge.writerow([relation[n][0], relation[n + 1][0], 'translated'])

            identifier = list(self.nodes.keys())
            name_and_metadata = list(self.nodes.values())
            for rel in self.rels:
                if query_from_sql:
                    from .models import query_data
                    if query_data(rel[0]):
                        value_list_1 = list(query_data(rel[0]).values())[0]
                        acc_num_1 = value_list_1[2]
                        tax_num_1 = value_list_1[3]
                        n_1 = name_and_metadata.index({rel[0]: 'Protein' + " " + f"{acc_num_1} /"+ f"Taxonomy:{tax_num_1}"})

                    if query_data(rel[2]):
                        value_list_2 = list(query_data(rel[2]).values())[0]
                        acc_num_2 = value_list_2[2]
                        tax_num_2 = value_list_2[3]
                        n_2 = name_and_metadata.index({rel[2]: 'Protein' + " " + f"{acc_num_2} /"+ f"Taxonomy:{tax_num_2}"})

                    tsv_writer_edge.writerow([identifier[n_1], identifier[n_2], rel[1]])
                else:
                    n_1 = name_and_metadata.index({rel[0]: 'Protein'+ " " + str(self.enrich_identifier_info[f"hgnc_{rel[0]}"]["UniProt ID"][0])})
                    n_2 = name_and_metadata.index({rel[2]: 'Protein'+ " " + str(self.enrich_identifier_info[f"hgnc_{rel[2]}"]["UniProt ID"][0])})
                    tsv_writer_edge.writerow([identifier[n_1], identifier[n_2], rel[1]])

    def enrich_generate_databaseinfo(self) -> tuple:
        hgnc_data = []
        uniprot_data = []
        id_hgnc = 1
        id_uniprot = 1
        for key, value in self.enrich_identifier_info.items():
            if key.startswith("hgnc"):
                hgnc_data.append((id_hgnc, value["HGNC ID"], value["Ensembl Gene ID"], value["UniProt ID"], key.split("_")[1]))
                id_hgnc += 1

            elif key.split("_")[0] == "uniprot":
                for i, y, z, uniprot, n in hgnc_data:
                    if key.split("_")[1] in uniprot:
                        uniprot_data.append((id_uniprot,
                                             value["first accession number"],
                                             value["NCBI taxonomy ID"],
                                             value["Name of protein"],
                                             value["full name of protein"],
                                             i))
                id_uniprot += 1
        new_hgnc_data = [(i, y, z, n) for i, y, z, uniprot, n in hgnc_data]
        return (new_hgnc_data, uniprot_data)

    def enrich_import_graph(self, edge_path: str, identifier: bool = False) -> nx.Graph:
        """
        Import the graph: nodes with HGNC symbol and their type (DNA, RNA and Protein)
        Parameters
        ----------
        edge_path: str
                  The .tsv file of enriched edge_list.
        Returns
        -------
        nx.Graph
        """
        Data = pd.read_csv(edge_path, sep='\t', names=['node_1', 'node_2', 'metadata'])

        # replace identifier with HGNC and type of molecule.
        mapping = {}
        if identifier:
            for rel in self.new_rels:
                if rel[2].split(" ")[0] == 'DNA':
                    mapping[int(rel[0])] = rel[2].split("/")[0]
                elif rel[2]== 'RNA':
                    mapping[int(rel[0])] = 'RNA'+ " "+ rel[1]
                elif rel[2].split(" ")[0] == 'Protein':
                    mapping[int(rel[0])] = rel[2]

        elif not identifier:
            for rel in self.new_rels:
                if len(rel[2]) > 3:
                    mapping[int(rel[0])] = rel[1] + " " + rel[2].split(" ")[0]
                else:
                    mapping[int(rel[0])] = rel[1] + " " + rel[2]
        new_data = Data.replace(mapping)
        logger.info("Import graph.")

        self.graph = nx.from_pandas_edgelist(new_data, source='node_1', target='node_2', edge_attr='metadata')


    def enrich_network(self, graph_output: str, print_edge_label: bool = False, identifier: bool = False) -> None:
        """
        Generate network of enriched info from imported enriched graph.
        Parameters
        ----------
        graph_output: str
                     The final output of a graph
        print_edge_label: bool
                     If true, will print the edge label of two nodes
        Returns
        -------
        None
        """
        self.check_output(graph_output)
        # specify colors for nodes
        color = {}
        if identifier:
            for node_id in self.graph.nodes():
                if node_id[:3] == 'DNA':
                    color[node_id] = 'green'
                elif node_id[:3] == 'RNA':
                    color[node_id] = 'blue'
                elif node_id[:7] == 'Protein':
                    color[node_id] = 'red'
        else:
            for node_id in self.graph.nodes():
                if node_id.split(" ")[1] == 'DNA':
                    color[node_id] = 'green'
                elif node_id.split(" ")[1] == 'RNA':
                    color[node_id] = 'blue'
                elif node_id.split(" ")[1] == 'Protein':
                    color[node_id] = 'red'

        node_colors = list(color.values())

        # specify color of edge
        edge_colors = 'black'

        # plot figure
        plt.figure(figsize=(18, 18))
        self.graph.pos = nx.spring_layout(self.graph, k=0.06)
        nx.draw_networkx(self.graph, pos=self.graph.pos,
                         node_size=150,
                         font_size=8,
                         with_labels=True,
                         node_color=node_colors,
                         edge_color=edge_colors,
                         alpha=0.3)
        if print_edge_label:
            edge_label = {(node1, node2): type['metadata'] for (node1, node2, type) in self.graph.edges(data=True)}
            nx.draw_networkx_edge_labels(self.graph, pos=self.graph.pos, edge_labels=edge_label)

        plt.savefig(graph_output)

    def add_data_to_database(self) -> None:
        from .models import add_data_hgnc, add_data_uniprot
        #a = Analyzer({}, None, ppi, None, None)
        relations = self.read_ppis(self.ppi_file)
        self.write_node_list("nodes_reduced.tsv")
        self.enrich_gather_identifier("nodes_reduced.tsv")
        hgnc_data = self.enrich_generate_databaseinfo()[0]
        uniprot_data = self.enrich_generate_databaseinfo()[1]
        add_data_hgnc(hgnc_data)
        add_data_uniprot(uniprot_data)
        logger.info("Database for HGNC and Uniprot are generated, and DATA is stored.")
        os.remove("nodes_reduced.tsv")


class Statistics(Network):
    def __init__(self,
                 nodes: dict,
                 graph: nx.Graph,
                 ppi_file: Optional[str],
                 node_path: Optional[str],
                 edge_path: Optional[str]):
        super().__init__(nodes, graph, ppi_file, node_path, edge_path)

    def enrich_generate_node_dict(self) -> dict:
        """
        A function used to get node dictionary of enriched info( DNA, RNA and Protein)
        Returns
        -------
        dict
        A dictionary of enriched info.
        """
        from collections import defaultdict
        self.nodes = defaultdict(dict)
        # A set contains unique node.
        nodes_ = set()
        for rel in self.rels:
            nodes_.add(rel[0])
            nodes_.add(rel[2])
        # A nested dictionary contains identifier as key, and in inner dict, the HGNC symbol as key.
        identifier = 1
        for symbol in nodes_:
            self.nodes[identifier][symbol] = 'DNA'
            identifier += 1
            self.nodes[identifier][symbol] = 'RNA'
            identifier += 1
            self.nodes[identifier][symbol] = 'Protein'
            identifier += 1

    def enrich_write_node_list(self, node_path: str) -> None:
        """
        Write node with their corresponding info( DNA, RNA, Protein) to the node_list file.
        Parameters
        ----------
        node_path: str
                    .tsv file with identifier and their info (DNA, RNA, Protein)
        Returns
        -------
        None
        """
        with open(node_path, 'wt') as outfile_1:
            tsv_writer_node = csv.writer(outfile_1, delimiter='\t')
            for identifier, info in self.nodes.items():
                for name, type in info.items():
                    tsv_writer_node.writerow([identifier, name, type])


    def enrich_node_label(self, relations: Iterable[Tuple[str]]) -> dict:
        """
        If there is node_list with only identifier and HGNC symbol,
        will generate node dictionary with other identifier and DNA, RNA and protein info.
        Parameters
        ----------
        relations: Iterable[Tuple[str]]
                  The relations from node list( only has identifer and HGNC info)
        Returns
        -------
        dict
        node dictionary with key as identifier and values as DNA/RNA/Protein info.
        """
        from collections import defaultdict
        self.nodes = defaultdict(dict)
        identifier = 1
        for index, symbol in relations:
            self.nodes[identifier][symbol] = 'DNA'
            identifier += 1
            self.nodes[identifier][symbol] = 'RNA'
            identifier += 1
            self.nodes[identifier][symbol] = 'Protein'
            identifier += 1



    def enrich_edge_from_ppi(self, relations: Iterable[Tuple[str]], edge_path: str) -> None:
        """
        Generate edge_list file if there is PPIs input
        Parameters
        ----------
        relations: Iterable[Tuple[str]]
                  The relationship from enriched node list file.
        edge_path: str
                  The .tsv file path of enriched edge list
        Returns
        -------
        None
        """

        with open(edge_path, 'wt') as outfile_2:
            tsv_writer_edge = csv.writer(outfile_2, delimiter='\t')
            for n in range(len(relations) - 1):
                if relations[n][2] == 'DNA':
                    tsv_writer_edge.writerow([relations[n][0], relations[n + 1][0], 'transcribed'])
                elif relations[n][2] == 'RNA':
                    tsv_writer_edge.writerow([relations[n][0], relations[n + 1][0], 'translated'])

            identifier = list(self.nodes.keys())
            name_and_metadata = list(self.nodes.values())
            for rel in self.rels:
                n_1 = name_and_metadata.index({rel[0]: 'Protein'})
                n_2 = name_and_metadata.index({rel[2]: 'Protein'})
                tsv_writer_edge.writerow([identifier[n_1], identifier[n_2], rel[1]])

    def open_original_edge(self, edge_path: str) -> Iterable[Tuple[str]]:
        """
        Function used to open original edge list with only identifier and HGNC symbol
        Parameters
        ----------
        edge_path: str
                  The .tsv file of edge list(identifier(protein) and type of interaction)
        Returns
        -------
        Iterable[Tuple[str]]
        """
        f = open(edge_path)
        content = f.readlines()
        self.edge_rels = [tuple(x.strip().split('\t')) for x in content]
        return self.edge_rels

    def enrich_edge_from_old_edge(self, relations: Iterable[Tuple[str]], new_edge_output: str) -> None:
        """
        Enrich edge list from original input edge list and
        Parameters
        ----------
        relations: Iterable[Tuple[str]]
                The relationship from enriched node list file
        new_edge_output: str
                The edge_list of enriched info.
        Returns
        -------
        None
        """
        with open(new_edge_output, 'wt') as outfile_2:
            tsv_writer_edge = csv.writer(outfile_2, delimiter='\t')
            for n in range(len(relations) - 1):
                if relations[n][2] == 'DNA':
                    tsv_writer_edge.writerow([relations[n][0], relations[n + 1][0], 'transcribed'])
                elif relations[n][2] == 'RNA':
                    tsv_writer_edge.writerow([relations[n][0], relations[n + 1][0], 'translated'])


            identifier = list(self.nodes.keys())
            name_and_metadata = list(self.nodes.values())
            for elem in self.edge_rels:
                n_1 = name_and_metadata.index({self.original_dict[elem[0]]:'Protein'})
                n_2 = name_and_metadata.index({self.original_dict[elem[1]]:'Protein'})
                tsv_writer_edge.writerow([identifier[n_1], identifier[n_2], elem[2]])


    def enrich_import_graph(self, edge_path: str) -> nx.Graph:
        """
        Import the graph: nodes with HGNC symbol and their type (DNA, RNA and Protein)
        Parameters
        ----------
        edge_path: str
                The .tsv file of enriched edge_list.
        Returns
        -------
        nx.Graph
        """
        Data = pd.read_csv(edge_path, sep='\t', names=['node_1', 'node_2', 'metadata'])

        mapping = {int(rel[0]): rel[1] + " " + rel[2] for rel in self.new_rels}
        new_data = Data.replace(mapping)

        self.graph = nx.from_pandas_edgelist(new_data, source='node_1', target='node_2', edge_attr='metadata')
        return self.graph

    def enrich_network(self, graph_output: str, print_edge_label: bool = True) -> None:
        """
        Generate network of enriched info from imported enriched graph.
        Parameters
        ----------
        graph_output: str
                    The final output of a graph
        print_edge_label: bool
                    If true, will print the edge label of two nodes
        Returns
        -------
        None
        """
        self.check_output(graph_output)
        # specify colors for nodes
        color = {}
        for node_id in self.graph.nodes():
            if node_id.split(" ")[1] == 'DNA':
                color[node_id] = 'green'
            elif node_id.split(" ")[1] == 'RNA':
                color[node_id] = 'blue'
            elif node_id.split(" ")[1] == 'Protein':
                color[node_id] = 'red'
        node_colors = list(color.values())

        # specify color of edge
        edge_colors = 'black'

        # plot figure
        plt.figure(figsize=(18, 18))
        self.graph.pos = nx.spring_layout(self.graph, k=0.06)
        nx.draw_networkx(self.graph, pos=self.graph.pos,
                         with_labels=True,
                         node_color=node_colors,
                         edge_color=edge_colors,
                         alpha=0.3)
        if print_edge_label:
            edge_label = {(node1, node2): type['metadata'] for (node1, node2, type) in self.graph.edges(data=True)}
            nx.draw_networkx_edge_labels(self.graph, pos=self.graph.pos, edge_labels=edge_label)

        plt.savefig(graph_output)


    def summary_statistics(self, enrich: bool = False) -> pandas.DataFrame:
        """
        Summary the statistics of network.
        Parameters
        ----------
        enrich: bool
                Control if it is enriched.
        Returns
        -------
        pandas.DataFrame
        The final summary of info in pandas file.
        """
        import math
        number_of_nodes = self.graph.number_of_nodes()
        number_of_edges = self.graph.number_of_edges()
        density = number_of_edges / (number_of_nodes - 1)
        local_node_connectivity = [nx.node_connectivity(self.graph, source, target) for source, target in self.graph.edges()]
        number_of_node_pair = math.factorial(number_of_nodes) / (math.factorial(number_of_nodes - 2) * 2)
        average_node_connectivity = float(sum(local_node_connectivity) / number_of_node_pair)
        if enrich:
            node_dna = []
            node_rna = []
            node_protein = []

            for node_id in self.graph.nodes():
                if node_id.split(" ")[1] == 'DNA':
                    node_dna.append(node_id)
                elif node_id.split(" ")[1] == 'RNA':
                    node_rna.append(node_id)
                elif node_id.split(" ")[1] == 'Protein':
                    node_protein.append(node_id)

            number_of_transcribed = 0
            number_of_translated = 0
            for node_1, node_2 in self.graph.edges():
                if node_1.split(" ")[1] == 'DNA' and node_2.split(" ")[1] == 'RNA':
                    number_of_transcribed += 1
                elif node_1.split(" ")[1] == 'RNA' and node_2.split(" ")[1] == 'Protein':
                    number_of_translated += 1
            number_of_ppis = number_of_edges - number_of_transcribed - number_of_translated

            d = {'Nodes': [number_of_nodes],
                 'Nodes-DNA': [len(node_dna)],
                 'Nodes-RNA': [len(node_rna)],
                 'Nodes-Protein': [len(node_protein)],
                 'Edges': [number_of_edges],
                 'Edges-transcribed': [number_of_transcribed],
                 'Edges-translated': [number_of_translated],
                 'Edges-PPI': [number_of_ppis],
                 'Density of Network': [density],
                 'Average node connectivity': [average_node_connectivity]}
        else:
            d = {'Nodes': [number_of_nodes],
                 'Nodes-DNA': 0,
                 'Nodes-RNA': 0,
                 'Nodes-Protein': [number_of_nodes],
                 'Edges': [number_of_edges],
                 'Edges-transcribed': 0,
                 'Edges-translated': 0,
                 'Edges-PPI': [number_of_edges],
                 'Density of Network': [density],
                 'Average node connectivity': [average_node_connectivity]}

        df = pd.DataFrame(data = d)
        return df


    def export_stats(self, data: pandas.DataFrame, stats_output: str) -> None:
        """
        Export the stats info to output file.
        Parameters
        ----------
        data: pandas.DataFrame
             Info from statistics function.
        stats_output: str
                     The output file.
        Returns
        -------
        None
        """
        export_basename = os.path.basename(stats_output)
        if export_basename.split(".")[1] == 'json':
            data.to_json("stats.json")
        elif export_basename.split(".")[1] == 'csv':
            data.to_csv("stats.csv", sep = '\t')
        elif export_basename.split(".")[1] == 'tsv':
            data.to_csv("stats.tsv", sep = '\t')
        elif export_basename.split(".")[1] == 'txt':
            data.to_csv('stats.txt', sep = '\t', mode = 'w')
        else:
            raise ValueError("Stats output path not csv, tsv, or json!")







