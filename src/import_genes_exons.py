""" import_genes_exons.py: parses gene and exon data from GFF3 formatted source and loads the records to a database. """

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
                        help='Path to the SQLite database to load')
    parser.add_argument('filepath',
                        help='Path to the gff3 file to load.')
    parser.add_argument('taxonid',
                        help='Taxon ID for the species being loaded.')
    parser.add_argument('-c', '--create', action='store_true',
                        help='Creates new database tables (after dropping any pre-existent tables).')
    args = parser.parse_args()
    return args


def create_tables(db_conn):
    """
    Creates "gene" and "exon" tables after dropping any pre-existing ones.

    :param db_conn: a connection to an sqlite3 database
    :return: None
    """
    cur = db_conn.cursor()
    cur.execute('''DROP TABLE IF EXISTS gene''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS gene (
            gene_id TEXT,
            gene_taxonid INTEGER,
            gene_symbol TEXT,
            gene_chr TEXT,
            gene_start_pos INTEGER,
            gene_end_pos INTEGER,
            gene_strand TEXT,
            gene_type TEXT,
            gene_name TEXT,
            PRIMARY KEY("gene_id","gene_taxonid","gene_chr")
        )
        ''')
    cur.execute('''CREATE INDEX gene_start_pos_idx ON gene(gene_taxonid, gene_chr, gene_start_pos)''')
    cur.execute('''CREATE INDEX gene_end_pos_idx ON gene(gene_taxonid, gene_chr, gene_end_pos)''')
    cur.execute('''CREATE INDEX gene_id_idx ON gene(gene_id)''')
    cur.execute('''CREATE INDEX gene_pos_idx ON gene(gene_chr, gene_start_pos, gene_end_pos)''')
    cur.execute('''CREATE INDEX gene_taxonid_symbol_idx ON gene(gene_taxonid, gene_symbol, gene_chr, gene_type)''')

    cur.execute('''DROP TABLE IF EXISTS exon''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS exon (
           exon_id TEXT,
           parent_gene TEXT,
           taxonid INTEGER,
           exon_chr TEXT,
           exon_start_pos INTEGER,
           exon_end_pos INTEGER,
           PRIMARY KEY("exon_id","taxonid","exon_chr")
        )
        ''')
    cur.execute('''CREATE INDEX exon_idx ON exon(parent_gene, exon_start_pos)''')
    db_conn.commit()


def load_genes_exons(db_conn, in_file, taxon_id):
    """
    Parse the GFF3 source and load the gene and exon data into the 'gene' and 'exon' tables.

    :param db_conn: a connection to an sqlite3 database
    :param in_file: path to the GFF3 source data file
    :param taxon_id: the species taxonomy ID for which genes/exons are loaded
    :return: None
    """
    gene_query = """INSERT INTO gene (
                     gene_id, gene_taxonid, gene_symbol, gene_chr, gene_start_pos, gene_end_pos, gene_strand, 
                     gene_type, gene_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                 """
    exon_query = """INSERT INTO exon (
                    exon_id, parent_gene, taxonid, exon_chr, exon_start_pos, exon_end_pos) VALUES 
                    (?,?,?,?,?,?)
                 """
    cur = db_conn.cursor()

    with flex_open(in_file) as in_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                # gene names are optional
                try:
                    gene_name = feature.qualifiers.get('Name')[0]
                except TypeError as e:
                    gene_name = None

                try:
                    gene_id = feature.qualifiers.get('Dbxref')[0]
                except TypeError as e:
                    logger.warning(f"\nwarning: GFF3 source data parsing failed - "
                                   f"'{feature.id}' is missing 'Dbxref' value")
                    continue

                rec_id = rec.id.replace("chr", "")

                gene_record = (gene_id,
                               taxon_id, feature.qualifiers.get('Symbol')[0], rec_id, feature.location.start.position,
                               feature.location.end.position, feature.location.strand, feature.type, gene_name)
                try:
                    cur.execute(gene_query, gene_record)
                except sqlite3.IntegrityError as e:
                    logger.error(f"\n'gene' table insertion failed: {gene_record}; \n{e}\n")
                    db_conn.rollback()
                    db_conn.close()
                    return False
                for exon in feature.sub_features:
                    exon_record = (exon.id, gene_id, taxon_id, rec_id, exon.location.start.position,
                                   exon.location.end.position)
                    try:
                        cur.execute(exon_query, exon_record)
                    except sqlite3.IntegrityError as e:
                        logger.error(f"\n'exon' table insertion failed: {exon_record}; \n{e}\n")
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

    status = load_genes_exons(db_conn, args.filepath, args.taxonid)
    if status is True:
        print(f"  Loading gene and exon records for <{args.taxonid}> into the database successfully completed..\n")
    else:
        print(f"  Loading gene and exon records for <{args.taxonid}> into the database failed - check the error log..\n")


if __name__ == '__main__':
    main()
