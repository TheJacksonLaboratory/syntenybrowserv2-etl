""" import_cytogenetic_loc.py: parses cytogenetic bands data from GFF3 source and loads the records to a database. """

import sys

import argparse
import logging
import sqlite3
from BCBio import GFF

from flex_open import flex_open

__version__ = "1.0.1"
__maintainer__ = "The Jackson Laboratory, Computational Sciences"
__email__ = "synbrowser-support@jax.org"
__status__ = "Production"


logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    # database name: required
    parser.add_argument('database',
                        help='Path to the database to load')
    parser.add_argument('filepath',
                        help='Path to the gff3 file to load.')
    parser.add_argument('taxonid',
                        help='Taxon ID for the cytogenetic locations being loaded.')
    parser.add_argument('-c', '--create', action='store_true',
                        help='Creates new database tables (after dropping any pre-existent tables).')
    args = parser.parse_args()
    return args


def create_tables(db_conn):
    """
    DB operations to drop and then create the 'cytogenetic_band' table.

    :param db_conn: a connection to an sqlite3 database
    :return: None
    """
    c = db_conn.cursor()
    c.execute('''DROP TABLE IF EXISTS cytogenetic_band''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS cytogenetic_band(
          id TEXT,
          taxon_id INTEGER,
          chr TEXT,
          source TEXT,
          type TEXT,
          start INTEGER,
          end INTEGER,
          location TEXT,
          color TEXT
        )
        ''')
    db_conn.commit()


def load_cytogenetic_band(db_conn, in_file, taxon_id):
    """
    Loads cytogenetic locations from a GFF3-format file: format is specified at http://gmod.org/wiki/GFF3.

    :param db_conn: a connection to an sqlite3 database
    :param in_file: the path to the file to be loaded
    :param taxon_id: the taxon ID for the locations we are loading
    :return: None
    """
    cyto_query = """INSERT INTO cytogenetic_band (
                       id, taxon_id, chr, source, type, start, end, 
                       location, color) 
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
                       """

    c = db_conn.cursor()
    with flex_open(in_file) as in_handle:
        for rec in GFF.parse(in_handle):
            chromosome = rec.id.replace('chr', '')
            for feature in rec.features:
                cyto_record = (feature.qualifiers.get('ID')[0], taxon_id, chromosome, feature.qualifiers.get('source')[0],
                               feature.type, feature.location.start.position+1, feature.location.end.position,
                               feature.qualifiers.get('Location')[0], feature.qualifiers.get('Color')[0])
                try:
                    c.execute(cyto_query, cyto_record)
                except sqlite3.IntegrityError as e:
                    logger.error(f"'cytogenetic_band' table insertion failed: {cyto_query}, {e}")
                    db_conn.rollback()
                    db_conn.close()
                    return False
    db_conn.commit()
    db_conn.close()
    return True


def main():
    args = parse_args()
    db_conn = sqlite3.connect(args.database)
    if args.create:
        create_tables(db_conn)

    status = load_cytogenetic_band(db_conn, args.filepath, args.taxonid)
    if status is True:
        print(f"  Loading band records for <{args.taxonid}> into the database successfully completed..\n")
    else:
        print(f"  Loading band records for <{args.taxonid}> into the database failed - check the error log..\n")


if __name__ == '__main__':
    main()
