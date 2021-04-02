""" import_homologs.py: parses homologs data from a tab separated file and loads the records into a database. """

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

# the source tab-separated file should have the following column names, where
# 'id' columns are the gene's official ID from the taxon authority, e.g. MGI for mouse;
# 'symbol' columns are the gene's symbols
HOM_FILE_HEADER_COLUMNS = [
    'type',
    'taxonid1', 'id1', 'symbol1', 'seqid1', 'start1', 'end1',
    'taxonid2', 'id2', 'symbol2', 'seqid2', 'start2', 'end2',
]


def parse_args():
    parser = argparse.ArgumentParser()
    # database name (required)
    parser.add_argument('database',
                        help='Path to the database to load')
    parser.add_argument('homologs',
                        help='Path to the tab separated source file to load.')
    parser.add_argument('-c', '--create', action='store_true',
                        help='Creates new database table (after dropping any pre-existent one).')
    args = parser.parse_args()
    return args


def create_table(db_conn):
    """
    Create the 'homolog' table after dropping ny pre-existing one.

    :param db_conn: a connection to an sqlite3 database
    :return: None
    """
    cur = db_conn.cursor()

    cur.execute('''DROP TABLE IF EXISTS homolog''')
    cur.execute('''
        CREATE TABLE homolog (
            ref_gene_id TEXT,
            ref_gene_sym TEXT,
            ref_taxon_id INTEGER,
            ref_seq_id TEXT,
            ref_start INTEGER,
            ref_end INTEGER,
            comp_gene_id TEXT,
            comp_gene_sym TEXT,
            comp_taxon_id INTEGER,
            comp_seq_id TEXT,
            comp_start INTEGER,
            comp_end INTEGER,
            PRIMARY KEY (ref_gene_id, ref_taxon_id, comp_gene_id, comp_taxon_id)
        )''')
    db_conn.commit()
    cur.execute('''CREATE INDEX homolog_comp_gene_id_idx ON
                    homolog(comp_gene_id, ref_gene_id)''')
    cur.execute('''CREATE INDEX homolog_ref_taxon_gene_idx ON
                    homolog(ref_taxon_id, ref_gene_id)''')
    cur.execute('''CREATE INDEX homolog_comp_taxon_gene_idx ON
                    homolog(comp_taxon_id, comp_gene_id)''')
    cur.execute('''CREATE INDEX homolog_ref_chr_idx ON
                    homolog(ref_seq_id)''')


def load_homologs(db_conn, homolog_filepath):
    """
    Load the contents of the file into the database.

    Each homolog is loaded twice, swapping reference and comparison organisms.
    :param db_conn: a connection to an sqlite3 database
    :param homolog_filepath: path to the file to be loaded
    :return num_records: number of homologs loaded (twice as many as rows in the file)
    """
    hom_file = flex_open(homolog_filepath)

    header = hom_file.readline()
    # Remove leading ##, if it exists.
    if header.startswith('##'):
        header = header[2:]
    columns = [x.strip().lower() for x in header.split('\t')]

    # check that we have all the expected headers
    missing = False
    for col in HOM_FILE_HEADER_COLUMNS:
        if col not in columns:
            logger.error(f"\nerror: expected header {col} is missing.\n")
            missing = True
    if missing:
        return None

    # let the user know if there are extra columns, but keep going
    for col in columns:
        if col not in HOM_FILE_HEADER_COLUMNS:
            logger.warning(f"\nwarning: ignoring extra column {col}.\n")

    cur = db_conn.cursor()
    reader = csv.DictReader(hom_file, fieldnames=HOM_FILE_HEADER_COLUMNS, delimiter='\t')

    # Now load all the rows.
    homolog_query = """INSERT OR REPLACE INTO homolog (
                 ref_gene_id, ref_gene_sym, ref_taxon_id, 
                 ref_seq_id, ref_start, ref_end, 
                 comp_gene_id, comp_gene_sym, comp_taxon_id,
                 comp_seq_id, comp_start, comp_end)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
    num_records = 0
    for row in reader:
        # remove "chr" from the start of the 'seqid' columns, if it is there
        if row['seqid1'].upper().startswith('CHR'):
            row['seqid1'] = row['seqid1'][3:]
        if row['seqid2'].upper().startswith('CHR'):
            row['seqid2'] = row['seqid2'][3:]

        if row['type'].strip().upper() != 'ORTHOLOGUE':
            logger.warning(f"\nwarning: unexpected type found. "
                           f"Expected 'orthologue', found {row['Type']}. "
                           f"Line is {row}.\n")

        # insert each relationship both ways
        reference = (row['id1'], row['symbol1'], row['taxonid1'],
                     row['seqid1'], row['start1'], row['end1'])
        comparison = (row['id2'], row['symbol2'], row['taxonid2'],
                      row['seqid2'], row['start2'], row['end2'])
        try:
            cur.execute(homolog_query, reference + comparison)
            cur.execute(homolog_query, comparison + reference)
            num_records = num_records + 2
        except sqlite3.IntegrityError as e:
            logger.error(f"\n'homolog' table insertion failed: {homolog_query}; \n{e}\n")
            db_conn.rollback()
            db_conn.close()
            return None
    db_conn.commit()
    db_conn.close()
    return num_records


def main():
    args = parse_args()
    db_conn = sqlite3.connect(args.database)
    if args.create:
        create_table(db_conn)

    num_records = load_homologs(db_conn, args.homologs)
    if num_records is not None:
        print(f"  Loading {num_records} homolog records into the database successfully completed..\n")
    else:
        print(f"  Loading homolog records into the database failed - check out the error log..\n")


if __name__ == '__main__':
    main()
