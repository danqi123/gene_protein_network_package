"This module is used to create tables for HGNC and UniProt information."

from .Utils import Profiler
import requests
import logging
from .startup import CONN_STRING
from sqlalchemy.orm import declarative_base, Session
from sqlalchemy import Column, Integer, String, create_engine, ForeignKey
import pymysql

logger = logging.getLogger('models')

engine = create_engine("mysql+pymysql://root:plab2rocks@localhost:3306/plab2db")
Base = declarative_base()


class Hgnc(Base):
    __tablename__ = "hgnc"

    id = Column(Integer, primary_key=True)
    HGNC_ID = Column(String(255))
    Ensembl_Gene_ID = Column(String(255))
    Gene_symbol = Column(String(255))


class Uniprot(Base):
    __tablename__ = "uniprot"

    id = Column(Integer, primary_key=True)
    accession_number = Column(String(255))
    NCBI_taxonomy_ID = Column(String(255))
    Name_of_protein = Column(String(255))
    full_name_protein = Column(String(255))
    HGNC_table_ID = Column(Integer, ForeignKey(Hgnc.id))

def delete_table():
    Base.metadata.drop_all(bind=engine)

def create_table():
    Base.metadata.create_all(bind=engine)

session = Session(bind = engine)


def add_data_hgnc(hgnc_data: list):
    """ Add data to hgnc table
    return: None
    """

    HGNC_to_add = []
    for hgnc in hgnc_data:
        ID, hgnc_ID, ensembl_Gene_ID, gene_Symbol = hgnc
        new_hgnc = Hgnc(id = ID, HGNC_ID = hgnc_ID, Ensembl_Gene_ID = ensembl_Gene_ID, Gene_symbol = gene_Symbol)
        # check if the input data is already in the table:
        if not bool(session.query(Hgnc).filter_by(id = ID).first()):
            HGNC_to_add.append(new_hgnc)

    session.add_all(HGNC_to_add)
    session.commit()

def add_data_uniprot(uniprot_data: list):
    """
    Add data to uniprot table
    return: None
    """
    UNIPROT_to_add = []
    for uniprot in uniprot_data:
        ID, acc_num, tax_id, name_HUMAN, full_name, table_id = uniprot
        new_uniprot = Uniprot(id=ID, accession_number = acc_num, NCBI_taxonomy_ID = tax_id,
                              Name_of_protein = name_HUMAN, full_name_protein = full_name,
                              HGNC_table_ID = table_id)
        # check if the data is already in the table:
        #if not bool(session.query(Uniprot).filter_by(id = ID).first()):
        UNIPROT_to_add.append(new_uniprot)

    session.add_all(UNIPROT_to_add)
    session.commit()

def query_data(gene_symbol: str) -> dict:
    """
    Query data from sql database.
    input: str
           gene symbol
    return: dict
    """
    from sqlalchemy import select
    from sqlalchemy.sql import exists

    if session.query(exists().where(Hgnc.Gene_symbol == gene_symbol)).scalar():
        filter_hgnc_id = select(Hgnc.id).filter_by(Gene_symbol = gene_symbol)
        hgnc_id = session.execute(filter_hgnc_id).all()
        filter_hgnc_ensembl = select(Hgnc.Ensembl_Gene_ID).filter_by(Gene_symbol = gene_symbol)
        hgnc_ensembl = session.execute(filter_hgnc_ensembl).all()

        logger.info(f"Gene Data of Gene Symbol {gene_symbol} is in Hgnc SQL database.")
        if session.query(exists().where(Uniprot.HGNC_table_ID == filter_hgnc_id)).scalar():
            filter_uniprot_acc = select(Uniprot.accession_number).filter_by(HGNC_table_ID=filter_hgnc_id)
            uniprot_acc = session.execute(filter_uniprot_acc).all()
            filter_uniprot_tax = select(Uniprot.NCBI_taxonomy_ID).filter_by(HGNC_table_ID=filter_hgnc_id)
            uniprot_tax = session.execute(filter_uniprot_tax).all()
            return_dict = {gene_symbol: [hgnc_id[0][0], hgnc_ensembl[0][0], uniprot_acc[0][0], uniprot_tax[0][0]]}
            logger.info(f"Protein Data of Gene Symbol {gene_symbol} is in Uniprot SQL database. ")
            return return_dict
        else:
            return_dict = {gene_symbol: [hgnc_id[0][0], hgnc_ensembl[0][0], None, None]}
            logger.error(f"No protein info of {gene_symbol} gene symbol.")
            return return_dict
    else:
        logger.error((f"No gene info of {gene_symbol} gene symbol."))
        return

if __name__ == "__main__":
    print(query_data("CDK1"))








