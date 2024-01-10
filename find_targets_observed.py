import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import sys
import argparse
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Postgres username, password, and database name
POSTGRES_ADDRESS = os.getenv("DB_HOST")  # Insert your DB address if it's not on Panoply
POSTGRES_PORT = os.getenv("DB_PORT")
POSTGRES_USERNAME = os.getenv("DB_USERNAME")
POSTGRES_PASSWORD = os.getenv("DB_PASSWORD")  # Change this to your Panoply/Postgres password
POSTGRES_DBNAME = 'trapum_web'  # Database name
# A long string that contains the necessary Postgres login information
postgres_str = ('mysql+pymysql://{username}:{password}@{ipaddress}:{port}/{dbname}'
.format(username=POSTGRES_USERNAME,
password=POSTGRES_PASSWORD,
ipaddress=POSTGRES_ADDRESS,
port=POSTGRES_PORT,
dbname=POSTGRES_DBNAME))
# Create the connection
cnx = create_engine(postgres_str)



file_type_pattern = "'MSGPS_S_%%'"



query_get_filterbank_files = '''
SELECT DISTINCT p.id AS pointing_id, t.source_name as target, p.utc_start, t.ra, t.dec 
FROM data_product dp
JOIN pointing p ON dp.pointing_id = p.id
JOIN target t ON p.target_id = t.id
WHERE t.source_name LIKE %s 
''' %file_type_pattern


printable_query = query_get_filterbank_files 
print(printable_query)

df = pd.read_sql_query(query_get_filterbank_files, con=cnx)

df.to_csv('mmgps_sband_observed.csv', index=False)


