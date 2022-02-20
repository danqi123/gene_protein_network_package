"""Test module constants."""
from pathlib import Path

THIS_DIR = Path(__file__).parent
THIS_PATH = THIS_DIR.joinpath("data")
ppi = THIS_PATH.joinpath('ppis.csv')
nodes = THIS_PATH.joinpath('node_list.tsv')
edges = THIS_PATH.joinpath('edge_list.tsv')
ppi_reduced = THIS_PATH.joinpath('ppis_reduced.csv')

nodes_enrich = THIS_PATH.joinpath("node_list_enrich.tsv")
edges_enrich = THIS_PATH.joinpath("edge_list_enrich.tsv")
graph = THIS_PATH.joinpath("graph.png")
#
# """Test module constants."""
# from pathlib import Path
#
# TEST_FOLDER = Path(__file__).parent
# TEST_DATA_DIR = TEST_FOLDER.joinpath("data")
# PPI_FILE_PATH = TEST_DATA_DIR.joinpath("ppis.csv")
# NODE_LIST_PATH = TEST_DATA_DIR.joinpath("node_list.tsv")
# EDGE_LIST_PATH = TEST_DATA_DIR.joinpath("edge_list.tsv")
# IMPORT_PPI_PATH = TEST_DATA_DIR.joinpath("import_ppi.csv")
