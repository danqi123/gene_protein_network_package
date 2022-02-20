"""Collection of tests for the network and analyzer."""
import pandas as pd
from plab2 import network
import os
import pytest
from .constants import ppi, nodes, edges, ppi_reduced, nodes_enrich, edges_enrich, graph



#
# ppi = "ppis.csv"
# ppi_reduced = "ppis_reduced.csv"
# nodes = "node_list.tsv"
# edges = "edge_list.tsv"

s_ppi = network.Statistics({}, None, ppi, None, None)
class TestNetwork:
    """Tests for the methods in parent class - Network."""
    def test_read_ppis(self):
        """Tests the read_ppis method."""
        n_ppi = network.Network({}, None, ppi, None, None)
        relation = n_ppi.read_ppis(ppi)
        assert isinstance(relation, list)
        assert ("USP14", "physical association", "AR") in relation

    def test_write_node_list_relations(self):
        """Tests the write_node_list and relations method."""
        n_ppi = network.Network({}, None, ppi, None, None)
        n_ppi.write_node_list(nodes)
        assert os.path.isfile(nodes)
        assert os.stat(nodes).st_size != 0
        rels = n_ppi.relations(nodes)
        assert isinstance(rels, list)

    def test_write_edge_list(self):
        """Tests the write_edge_list method."""
        n_ppi = network.Network({}, None, ppi, None, None)
        n_ppi.write_edge_list(edges)
        assert os.path.isfile(edges)
        assert os.stat(edges).st_size != 0

    def test_check_output(self):
        """Tests the check_output methods."""
        n_ppi = network.Network({}, None, ppi, None, None)
        output = "graph.txt"
        with pytest.raises(ValueError):
            n_ppi.check_output(output)

    def test_import_graph(self):
        """Tests the function import graph."""
        # read ppis
        n_ppi = network.Network({}, None, ppi, None, None)
        number_of_edge = len(n_ppi.read_ppis(ppi)[1:])
        n_ppi.write_node_list(nodes)
        n_ppi.write_edge_list(edges)
        graph_edge_list = [edge for edge in n_ppi.import_graph(edges).edges]
        result_of_edges = len(graph_edge_list)
        assert number_of_edge == result_of_edges

        # read node and edge lists
        n_node_edge = network.Network({}, None, None, nodes, edges)
        n_node_edge.relations(nodes)
        number_of_edge_ = len(n_node_edge.relations(edges)[1:])
        graph_edge_list_ = [edge for edge in n_node_edge.import_graph(edges).edges]
        result_of_edges_ = len(graph_edge_list_)
        assert number_of_edge_ == result_of_edges_

    def test_generate_image(self):
        """Tests the function generate_graph_network. """
        if os.path.isfile(graph):
            os.remove(graph)
        assert not os.path.isfile(graph)

        n_ppi = network.Network({}, None, ppi, None, None)
        n_ppi.import_graph(edges)
        n_ppi.generate_graph_network(graph)
        assert os.path.isfile(graph)


class TestAnalyzer:
    """Test for functions in Analyzer class."""

    def test_shortest_path(self):
        """Tests the function shortest path."""
        source = "CREBBP"
        target = "TRA2B"
        path = ['CREBBP', 'CREB1', 'EP300', 'REL', 'BTRC', 'TP63', 'SUMO1',
                 'MTOR', 'GNB1', 'PIK3R1', 'BRCA1', 'BARD1', 'PGAM5', 'PRKN',
                 'PSMD1', 'BLM', 'RPA1', 'PAN2', 'SRSF3', 'S100A9', 'TRA2B']
        a_ppi = network.Analyzer({}, None, ppi, None, None)
        a_ppi.write_node_list(nodes)
        a_ppi.write_edge_list(edges)
        a_ppi.import_graph(edges)
        result_path = a_ppi.shortest_path(source, target)
        assert result_path is not None
        assert path in result_path

#     def test_import_graph(self):
#         """Tests the function import graph."""
#         pass
#
#
class TestStatistics:
    """Tests for function if Statistics class."""
    def test_enrich_graph(self):
        """Tests enrich the graph with RNA/DNA nodes and the associated edges."""
        s_ppi = network.Statistics({}, None, ppi, None, None)
        s_ppi.enrich_generate_node_dict()
        s_ppi.enrich_write_node_list(nodes_enrich)
        s_ppi.enrich_edge_from_ppi(s_ppi.relations(nodes_enrich), edges_enrich)
        graph_ = s_ppi.enrich_import_graph(edges_enrich)

        node_dna = []
        node_rna = []
        node_protein = []

        for node_id in graph_.nodes():
            if node_id.split(" ")[1] == 'DNA':
                node_dna.append(node_id)
            elif node_id.split(" ")[1] == 'RNA':
                node_rna.append(node_id)
            elif node_id.split(" ")[1] == 'Protein':
                node_protein.append(node_id)

        number_of_transcribed = 0
        number_of_translated = 0
        for node_1, node_2 in graph_.edges():
            if node_1.split(" ")[1] == 'DNA' and node_2.split(" ")[1] == 'RNA':
                number_of_transcribed += 1
            elif node_1.split(" ")[1] == 'RNA' and node_2.split(" ")[1] == 'Protein':
                number_of_translated += 1
        assert len(node_dna) == len(node_rna) == len(node_protein) == number_of_translated == number_of_transcribed

    def test_summary_statistics(self):
        """Tests function summary_statistics: if it can calculate summary statistics for a graph."""
        s_ppi = network.Statistics({}, None, ppi, None, None)
        s_ppi.enrich_generate_node_dict()
        s_ppi.enrich_write_node_list(nodes_enrich)
        s_ppi.enrich_edge_from_ppi(s_ppi.relations(nodes_enrich), edges_enrich)
        s_ppi.enrich_import_graph(edges_enrich)
        d = s_ppi.summary_statistics(enrich=True)
        assert isinstance(d, pd.DataFrame)

        nodes_dna = d.iloc[0, 1]
        nodes_rna = d.iloc[0, 2]
        nodes_protein = d.iloc[0, 3]
        transcribed = d.iloc[0, 5]
        translated = d.iloc[0, 6]
        assert nodes_dna == nodes_rna == nodes_protein == transcribed == translated == 621

    def test_export_stats(self):
        """Tests function export_stats."""
        s_ppi = network.Statistics({}, None, ppi, None, None)
        s_ppi.enrich_generate_node_dict()
        s_ppi.enrich_write_node_list(nodes_enrich)
        s_ppi.enrich_edge_from_ppi(s_ppi.relations(nodes_enrich), edges_enrich)
        s_ppi.enrich_import_graph(edges_enrich)
        d = s_ppi.summary_statistics(enrich=True)
        s_ppi.export_stats(d, "stats.json")
        assert os.path.isfile("stats.json")
        assert os.stat("stats.json").st_size != 0
        os.remove("stats.json")
        assert not os.path.isfile("stats.json")

        s_ppi.export_stats(d, "stats.csv")
        assert os.path.isfile("stats.csv")
        assert os.stat("stats.csv").st_size != 0
        os.remove("stats.csv")
        assert not os.path.isfile("stats.csv")

        s_ppi.export_stats(d, "stats.tsv")
        assert os.path.isfile("stats.tsv")
        assert os.stat("stats.tsv").st_size != 0
        os.remove("stats.tsv")
        assert not os.path.isfile("stats.tsv")

        with pytest.raises(ValueError):
            s_ppi.export_stats(d, "stats.png")



