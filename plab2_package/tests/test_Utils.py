"""Tests for the Utils."""

import os
from plab2.Utils import Profiler
from plab2.startup import HGNC_data_path, UniProt_data_path
HGNC_symbol = "RPL10"
accession_number = "O00483"
p = Profiler(HGNC_symbol)

class TestUtils:
    """Test for functions related to download and cache files."""
    def test_request(self):
        """Tests function for downloading data from HGNC."""
        p.request()
        path_root = HGNC_data_path
        file_name = f"/{HGNC_symbol}.json"
        path = path_root + file_name
        assert os.path.isfile(path)
        assert os.stat(path).st_size != 0

    def test_request_uniprot(self):
        """Tests function for downloading data from UniProt."""
        p.request_uniprot(accession_number)
        path_root_ = UniProt_data_path
        file_name_ = f"/{accession_number}.fasta"
        path_ = path_root_ + file_name_
        assert os.path.isfile(path_)
        assert os.stat(path_).st_size != 0

    def test_extract(self):
        """Test for extract HGNC info from cache file."""
        extract_info = p.extract()
        assert extract_info["HGNC ID"] == "HGNC:10298"
        assert extract_info["Ensembl Gene ID"] == "ENSG00000147403"
        assert extract_info["UniProt ID"] == ["P27635"]

    def test_extract_uniprot(self):
        """Test for extract UniProt info from cache file."""
        extract_info_ = p.extract_uniprot(accession_number)
        assert extract_info_["Name of protein"] == "NDUA4_HUMAN"
        assert extract_info_["full name of protein"] == "Cytochrome c oxidase subunit NDUFA4"
        assert extract_info_["NCBI taxonomy ID"] == "9606"
















