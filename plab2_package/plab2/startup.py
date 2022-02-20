import os
from pathlib import Path
import logging
from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy_utils.functions import database_exists


# identify home directory
home_dir = Path.home()


# join paths
HGNC_data_path = os.path.join(str(home_dir), ".wangd0", "data", "HGNC")
UniProt_data_path = os.path.join(str(home_dir), ".wangd0", "data", "UniProt")
logs_path = os.path.join(str(home_dir), ".wangd0", "logs")


# create folder
os.makedirs(logs_path, mode=0o777, exist_ok=True)
os.makedirs(HGNC_data_path, mode=0o777, exist_ok=True)
os.makedirs(UniProt_data_path, mode=0o777, exist_ok=True)

# create SQLite database
database_path = os.path.join(str(home_dir), ".wangd0", "plab2.db")
CONN_STRING = "sqlite:///" + database_path
if not database_exists(CONN_STRING):
    engine = create_engine(CONN_STRING)
    Base = declarative_base()
    Base.metadata.create_all(bind=engine)



# if use pathlib:
# ssh_path = f"{os.getenv('HOME')}/temp/.ssh")
# ssh = Path(ssh_path)
# ssh.mkdir(parents=true)

filepath = os.path.join(logs_path, "log_file_name.log")
logging.basicConfig(filename = filepath,
                    level = logging.INFO,
                    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')




