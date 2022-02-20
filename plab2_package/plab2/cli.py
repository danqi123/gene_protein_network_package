import click
import os
from .network import Network, Analyzer, Statistics
import logging
from .startup import HGNC_data_path, UniProt_data_path, CONN_STRING


logger = logging.getLogger('cli')


@click.group(help = f'The Command Line Utilities of generating network.')
def main():
    pass

@main.command()
@click.option('-s', '--source')
@click.option('-t', '--target')
@click.argument('output_path')
@click.option('-p', '--ppi', default=None, help="A CSV file containing PPIs.")
@click.option('-n', '--nodes', default=None, help="A TSV file containing defined nodes of a network.")
@click.option('-e', '--edges', default=None, help="A TSV file containing defined edges of a network.")
@click.option('-v', '--verbose', default=False, is_flag=True, help="When used, will print the paths to STDOUT.")
@click.option('--add_edge', default=False, is_flag=True)
def path(output_path:str, source:str, target:str, ppi: str, nodes:str, edges:str, verbose:bool, add_edge:bool):
    if ppi:
        logger.info("PPI file is accepted as input.")
        node_path, edge_path = "node_list.tsv", "edge_list.tsv"
    else:
        logger.info("node/edge lists are accepted as input.")
        node_path, edge_path = nodes, edges

    a = Analyzer({}, None, ppi, None, None)
    a.write_node_list(node_path)
    a.write_edge_list(edge_path)
    a.import_graph(edge_path)
    try:
        a.shortest_path(source, target, print_option = verbose)
    except:
        logger.warning(f'No path found between {source} and {target}!')
    a.generate_graph_network(output_path, add_edge)
    logger.info(f"New graph image was generated and its location is {output_path}")


@main.command()
@click.argument('ppi')
@click.argument('node_file')
@click.argument('edge_file')
@click.option('--enrich', default = False, is_flag = True, help = 'when used, will enrich the network')
@click.option('--query_from_sql', default = False, is_flag = True, help = 'when used, will query from database.')
def compile(ppi: str, node_file: str, edge_file: str, enrich: bool, query_from_sql: bool) -> None:
    # No enrichment of network
    if not enrich:
        n = Network({}, None, ppi, None, None)
        n.write_node_list(node_file)
        n.write_edge_list(edge_file)
        logger.info(f"New node/edge files were made and their locations are {node_file} and {edge_file}")

    # enrich network with gene and protein info.
    elif enrich:
        a = Analyzer({}, None, ppi, None, None)
        if query_from_sql:
            #import models
            a.enrich_generate_node_dict(query_from_sql)
            a.enrich_write_node_list(node_file)
            a.enrich_write_edge_list(a.relations(node_file), edge_file, query_from_sql)
        else:
            a.write_node_list("nodes_reduced.tsv")
            a.enrich_gather_identifier("nodes_reduced.tsv")
            a.enrich_generate_node_dict()
            a.enrich_write_node_list(node_file)
            a.enrich_write_edge_list(a.relations(node_file), edge_file)
            os.remove("nodes_reduced.tsv")
        logger.info(f"New node/edge files were made and their locations are {node_file} and {edge_file}")


@main.command()
@click.argument('ppi')
@click.argument('output')
@click.argument('nodes')
@click.argument('edges')
@click.option('-v', '--verbose', default = False, is_flag = True,
              help = 'when used, will print the edge label.')
@click.option('--enrich', default = False, is_flag = True,
              help = 'when used, will enrich the network')
@click.option('-s', '--show_identifier', default = False, is_flag = True,
              help = 'when used, will show identifier instead of HGNC symbol.')
@click.option('--query_from_sql', default = False, is_flag = True,
              help = 'when used, will query from database.')
def create(ppi: str, nodes: str, edges: str, output: str, verbose: bool, enrich: bool, show_identifier: bool, query_from_sql: bool):
    if not enrich:
        n = Network({}, None, ppi, None, None)
        relations = n.read_ppis(ppi)
        n.write_node_list(nodes)
        n.write_edge_list(edges)
        logger.info("New node/edge files were made and their locations are '/Exercise_5/node_list.tsv' and '/Exercise_5/edge_list.tsv'.")
        n.import_graph(edges)
        n.generate_graph_network(output, verbose)
        logger.info("Graph image was generated and its location is '/Exercise_5/graph.png'.")
    elif enrich:
        a = Analyzer({}, None, ppi, None, None)
        if query_from_sql:

            a.enrich_generate_node_dict(query_from_sql)
            a.enrich_write_node_list(nodes)
            a.enrich_write_edge_list(a.relations(nodes), edges, query_from_sql)
        else:
            relations = a.read_ppis(ppi)
            a.write_node_list("nodes_reduced.tsv")
            a.enrich_gather_identifier("nodes_reduced.tsv")
            a.enrich_generate_node_dict()
            a.enrich_write_node_list(nodes)
            a.enrich_write_edge_list(a.relations(nodes), edges)
            os.remove("nodes_reduced.tsv")
        a.enrich_import_graph(edges, identifier=show_identifier)
        a.enrich_network(output, verbose, identifier=show_identifier)
        logger.info("network which is shown in graph.")

@main.command()
@click.option('-p', '--ppi', default = None)
@click.option('-n','--node_file', default = None)
@click.option('-e','--edge_file', default = None)
@click.option('--export', default = None, help='when used, EXPORT the statistics to file.')
@click.option('--print_table', default = False, is_flag = True, help='When used, print table in STDOUT.')
@click.option('--enrich', default = False, is_flag = True, help='when used, enrich the DNA and RNA info in network.')


def stats(ppi: str, node_file: str, edge_file: str, enrich: bool, print_table: bool, export: str):
    if ppi and not node_file and not edge_file:
        logger.info("PPI file is accepted as input.")
        s = Statistics({}, None, ppi, None, None)
        if enrich:
            s.enrich_generate_node_dict()
            s.enrich_write_node_list('node_list_enrich.tsv')
            logger.info("New enriched node file was made and the location is 'node_list_enrich.tsv.")
            s.enrich_edge_from_ppi(s.relations('node_list_enrich.tsv'), 'edge_list_enrich.tsv')
            logger.info("New enriched edge file was made and the location is 'edge_list_enrich.tsv.")
            s.enrich_import_graph('edge_list_enrich.tsv')
            data = s.summary_statistics(enrich)
        elif not enrich:
            s.write_node_list('node_list.tsv')
            logger.info("New node file was made and the location is 'node_list.tsv.")
            s.write_edge_list('edge_list.tsv')
            logger.info("New edge file was made and the location is 'edge_list.tsv.")
            s.import_graph('edge_list.tsv')
            data = s.summary_statistics(enrich)
        if print_table:
            from tabulate import tabulate
            click.echo(tabulate(data, headers = 'keys', tablefmt = 'psql'))
        if export:
            s.export_stats(data, export)

    elif not ppi and node_file and edge_file:
        logger.info("node/edge lists are accepted as input.")
        s = Statistics({}, None, None, node_file, edge_file)
        relations = s.relations(node_file)
        if enrich:
            s.enrich_node_label(relations)
            s.enrich_write_node_list('node_list_enrich.tsv')
            logger.info("New enriched node file was made and the location is '/Exercise_5/node_list_enrich.tsv.")
            s.open_original_edge(edge_file)
            s.enrich_edge_from_old_edge(s.relations('node_list_enrich.tsv'),'edge_list_enrich.tsv')
            logger.info("New enriched edge file was made and the location is '/Exercise_5/edge_list_enrich.tsv.")
            s.enrich_import_graph('edge_list_enrich.tsv')
            data = s.summary_statistics(enrich)
        elif not enrich:
            s.import_graph(edge_file)
            data = s.summary_statistics(enrich)
        if print_table:
            from tabulate import tabulate
            click.echo(tabulate(data, headers='keys', tablefmt='psql'))
        if export:
            s.export_stats(data, export)

@main.command()
@click.argument("hgnc_symbol")
@click.option('--query_from_sql', default = False, is_flag = True, help='When used, query data from SQL.')
def info(hgnc_symbol: str, query_from_sql: bool = False):
    HGNC_root = "http://rest.genenames.org/fetch/symbol/"
    UNIPROT_root = "https://www.uniprot.org/uniprot/"
    if query_from_sql:
        from sqlalchemy.orm import declarative_base, Session
        from sqlalchemy import MetaData, Table, Column, Integer, String, create_engine, ForeignKey
        from .models import query_data
        engine = create_engine(CONN_STRING, echo=True)
        session = Session(bind=engine)
        if query_data(hgnc_symbol):
            list_info = list(query_data(hgnc_symbol).values())
            click.echo(f"HGNC ID is {list_info[0][0]}")
            click.echo(f"Ensembl Gene ID is {list_info[0][1]}")
            click.echo(f"Gene Symbol is {hgnc_symbol}")

            HGNC_id = f"{hgnc_symbol}"
            query = HGNC_root + HGNC_id
            click.echo(f"The HGNC link to this HGNC symbol: {query}")
            click.echo(f"UniProt accession number is {list_info[0][2]}")
            click.echo(f"NCBI Taxonomy ID is {list_info[0][3]}")

            UNIPROT_id = f"{list_info[0][2]}"
            query_uniprot = UNIPROT_root+UNIPROT_id
            click.echo(f"The Uniprot link to this HGNC symbol: {query_uniprot}")
    else:
        from .Utils import Profiler
        p = Profiler(hgnc_symbol)
        file_path = os.path.join(HGNC_data_path, f'{hgnc_symbol}.json')
        if not os.path.isfile(file_path):
            p.request()
        click.echo(f"The identifiers of {hgnc_symbol} are:")
        click.echo(p.extract())
        click.echo(f"The HGNC link to this HGNC symbol: {HGNC_root+hgnc_symbol}")
        for acc_num in p.extract()["UniProt ID"]:
            if acc_num is not None:
                uniprot_path_root = UniProt_data_path
                uniprot_file_name = f"/{acc_num}.fasta"
                uniprot_path = uniprot_path_root + uniprot_file_name
                if not os.path.isfile(uniprot_path):
                    p.request_uniprot(acc_num)
                click.echo(p.extract_uniprot(acc_num))
                click.echo(f"The UniProt link to this HGNC symbol: {UNIPROT_root+acc_num}")


@main.command()
@click.option('-p', '--ppi', default = None)
@click.option('--enrich', default = False, is_flag = True, help="when used, will enrich HGNC, ensembl, and UniProt info.")

def populate(ppi: str, enrich: bool):
    if ppi:
        logger.info("PPI file is accepted as input.")
        if enrich:
            #from .models import add_data_hgnc, add_data_uniprot
            a = Analyzer({}, None, ppi, None, None)
            a.add_data_to_database()


if __name__ == '__main__':
    main()
