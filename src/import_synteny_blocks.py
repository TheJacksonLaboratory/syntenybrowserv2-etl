""" import_synteny_blocks.py: parses syntenic blocks data from custom formatted source and loads the records into a DB. """

import sys

import argparse
import csv
import logging
import sqlite3

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
                        help='Path to the SQLite database to load')
    parser.add_argument('filepath',
                        help='Path to the tab delimited file to load.')
    parser.add_argument('-c', '--create', action='store_true',
                        help='Creates new database tables (after dropping any pre-existent tables).')

    args = parser.parse_args()
    return args


def create_tables(db_conn):
    """
    DB operation to drop and then create the 'syntenic_block' table.

    :param db_conn: a connection to an sqlite3 database
    :return: None
    """
    c = db_conn.cursor()
    c.execute('''DROP TABLE IF EXISTS syntenic_block''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS syntenic_block (
            ref_taxonid INTEGER,
            ref_chr TEXT,
            ref_start_pos INTEGER,
            ref_end_pos INTEGER,
            comp_taxonid INTEGER,
            comp_chr TEXT,
            comp_start_pos INTEGER,
            comp_end_pos INTEGER,
            same_orientation BOOLEAN,
            symbol TEXT,
            PRIMARY KEY (ref_taxonid, comp_taxonid, ref_chr, ref_start_pos))
    ''')

    c.execute('''CREATE UNIQUE INDEX IF NOT EXISTS syntenic_taxons_ref_idx ON
                      syntenic_block (ref_taxonid, comp_taxonid, ref_chr)''')
    db_conn.commit()


def validate_row_structure(record):
    """
    Runs basic validation on each record(row) from the source data.
    """
    assert(len(record) == 10), 'invalid data items number in input file'

    block_id = record[9]
    assert(block_id[0:11] == 'ID=SynBlock'), 'invalid block id format in input file'


def load_syntenic_blocks(db_conn, in_file):
    """
    Parse the source file and load the syntenic blocks data into the 'syntenic_block' table.

    :param db_conn: a connection to an sqlite3 database
    :param in_file: path to the source data file
    :return: None
    """
    c = db_conn.cursor()

    block_query = """INSERT OR REPLACE INTO syntenic_block (
                         ref_taxonid, ref_chr, ref_start_pos, ref_end_pos,
                         comp_taxonid, comp_chr, comp_start_pos, comp_end_pos,
                         same_orientation, symbol
                         )
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""

    with flex_open(in_file) as in_handle:
        csv_reader = csv.reader(in_handle, delimiter='\t')

        for row in csv_reader:
            if row:
                same_orientation = row[8] == '+'

                try:
                    validate_row_structure(row)

                    record_forward = (row[1], row[0], row[2], row[3], row[5], row[4],
                                      row[6], row[7], same_orientation, row[9][3:])
                    c.execute(block_query, record_forward)
                except AssertionError as ae:
                    logger.error(f"\nerror: data validation failed: {row}; \n{ae}\n")
                    return False
                except sqlite3.IntegrityError as e:
                    logger.error(f"\n'syntenic_block' table insertion failed: {record_forward}; \n{e}\n")
                    db_conn.rollback()
                    db_conn.close()
                    return False

                try:
                    record_reverse = (row[5], row[4], row[6], row[7], row[1], row[0],
                                      row[2], row[3], same_orientation, row[9][3:])
                    c.execute(block_query, record_reverse)
                except sqlite3.IntegrityError as e:
                    logger.error(f"\n'syntenic_block' table insertion failed: {record_forward}; \n{e}\n")
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

    status = load_syntenic_blocks(db_conn, args.filepath)
    if status is True:
        print(f"  Loading syntenic blocks into the database successfully completed..\n")
    else:
        print(f"  Loading syntenic blocks into the database failed - check the error log..\n")


if __name__ == '__main__':
    main()
