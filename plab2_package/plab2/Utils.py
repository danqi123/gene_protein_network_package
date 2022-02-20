"""Methods related to downloading and storing information HGNC and UniProt."""

import requests
import logging
from .startup import HGNC_data_path, UniProt_data_path

logger = logging.getLogger('Utils')

class Profiler():
    def __init__(self, protein_id: str, get_uniprot: bool = False):
        self.protein_id = protein_id
        self.get_uniprot = get_uniprot

    def request(self) -> None:
        """
        Use HGNC symbol as input, and will write corresponding .xml files using API method to data folder.
        ------
        Parameter: protein_id : str
                  HGNC symbol
        ------
        Return: None
               .xml file
        """
        headers = {"Accept": "application/json", }
        HGNC_root = "http://rest.genenames.org/fetch/symbol/"
        HGNC_id = f"{self.protein_id}"
        query = HGNC_root + HGNC_id
        r = requests.get(query, headers=headers)
        if r.ok:
            logger.info(f"Successful request for {self.protein_id}.")


            path_root = HGNC_data_path
            file_name = f"/{self.protein_id}.json"
            path = path_root + file_name
            content = r.json()["response"]
            if content["numFound"] > 0:
                with open(path, "w") as file:
                    file.write(r.text)
                logger.info(f"{self.protein_id}.json file is stored in data folder.")
        else:
            logger.warning(f"no results for {self.protein_id} were found")

    # # the status code is unacceptable.
    # else:
    #     logger.error(f"Request error: {protein_id} is not acceptable.")

    def extract(self) -> dict:
        """
            Extract HGNC ID, ENsembl Gene ID and Uniprot ID info to a dictionary.
            ------
            Parameter: prtein_id: str
                      HGNC symbol in PPI file
            ------
            Return: dict
                   return a dictionary contains info of HGNC id, ensembl id and uniprot id.
                   Also prints out the link to this HGNC symbol.
            """
        import json
        input_file = HGNC_data_path + f"/{self.protein_id}.json"

        file = open(input_file)
        f = file.read()
        data = json.loads(f)

        info_dict = {}

        docs = data["response"]["docs"][0]
        try:
            info_dict["HGNC ID"] = docs["hgnc_id"]
        except:
            info_dict["HGNC ID"] = [None]
        try:
            info_dict["Ensembl Gene ID"] = docs["ensembl_gene_id"]
        except:
            info_dict["Ensembl Gene ID"] = [None]
        try:
            info_dict["UniProt ID"] = docs["uniprot_ids"]
        except:
            info_dict["UniProt ID"] = [None]

        logger.info(f"Identifier info of {self.protein_id} is generated.")

        return info_dict

    def request_uniprot(self, accession_number: str) -> None:
        """ request info from UniProt.
        Input: accession number.
        """
        uniprot_root = "https://www.uniprot.org/uniprot/"
        full_url = uniprot_root + accession_number + ".fasta"
        r = requests.get(full_url)
        if r.ok:
            logger.info(f"Successful Uniprot request for {accession_number}.")
            path_root = UniProt_data_path
            file_name = f"/{accession_number}.fasta"
            path = path_root + file_name
            with open(path, "w") as file:
                file.write(r.text)
            logger.info(f"{self.protein_id}.fasta file is stored in data folder.")
        else:
            logger.warning("something went wrong", r.status_code)

    def extract_uniprot(self, accession_number: str) -> dict:
        """
        extract uniprot info from fasta file in cache.
        Input: accession number.
        """
        fasta_file = UniProt_data_path + f"/{accession_number}.fasta"
        f = open(fasta_file)
        content = f.readlines()
        info_string = content[0].lstrip(">sp|")
        first_acc_num = info_string.split("|")[0]
        other_info = info_string.split("|")[1]
        rest_info = other_info.split(" ")
        name_of_protein = rest_info[0]
        for elem in rest_info:
            if elem.startswith("OX="):
                NCBI_ID = elem[3:]
            elif elem.startswith("GN="):
                primary_gene_symbol = elem[3:]
            elif elem.startswith("OS="):
                OS_index = rest_info.index(elem)
                full_name_protein = " ".join(rest_info[1:OS_index])
        protein_info = {"first accession number": first_acc_num,
                        "primary gene symbol": primary_gene_symbol,
                        "NCBI taxonomy ID": NCBI_ID,
                        "Name of protein": name_of_protein,
                        "full name of protein": full_name_protein}
        logger.info(f"UniProt info of {accession_number} is generated.")
        return protein_info

#
#
# if __name__ == "__main__":
#     # p = Profiler("ATF2")
#     # p.request()
#     # for acc_num in p.extract()["UniProt ID"]:
#     #     p.request_uniprot(acc_num)
#     #     print(p.extract_uniprot(acc_num))
