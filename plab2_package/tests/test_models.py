"""Test for models."""

import os
import sqlite3
from plab2 import models
from plab2.startup import database_path
from plab2.startup import CONN_STRING
from sqlalchemy.orm import Session
from sqlalchemy import Column, Integer, String, create_engine, ForeignKey

engine = create_engine(CONN_STRING, echo = True)
session = Session(bind=engine)

class Testmodel:
    def test_create_database(self):
        """Tests whether a new database file is generated in the project directory."""
        assert os.path.isfile(database_path)

    def test_correct_schema_hgnc(self):
        """Tests the correct schema is generated."""
        conn = sqlite3.connect(database_path)
        c = conn.cursor()
        # get the count of tables with the name
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='hgnc' ''')
        # if the count is 1, then table exists
        assert c.fetchone()[0] == 1
        # check column exists
        assert session.query(models.Hgnc.id, models.Hgnc.Ensembl_Gene_ID, models.Hgnc.Gene_symbol, models.Hgnc.HGNC_ID)

    def test_correct_schema_uniprot(self):
        """Tests the correct schema is generated."""
        conn = sqlite3.connect(database_path)
        c = conn.cursor()
        # get the count of tables with the name
        c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='uniprot' ''')
        # if the count is 1, then table exists
        assert c.fetchone()[0] == 1
        # check column exists
        assert session.query(models.Uniprot.id, models.Uniprot.accession_number, models.Uniprot.NCBI_taxonomy_ID,
                             models.Uniprot.full_name_protein, models.Uniprot.HGNC_table_ID, models.Uniprot.Name_of_protein)

    def test_number_of_row(self):
        """Tests the table contains the expected number of rows. """
        # HGNC table
        query = session.query(models.Hgnc)
        nrows = query.count()
        assert nrows == 620
        # UniProt table
        query_uniprot = session.query(models.Uniprot)
        nrows_uniprot = query_uniprot.count()
        assert nrows_uniprot == 606

    def test_query_data(self):
        """Tests the function query data."""
        result = models.query_data("CDK1")
        assert result == {'CDK1': [66, 'ENSG00000170312', 'P06493', '9606']}
