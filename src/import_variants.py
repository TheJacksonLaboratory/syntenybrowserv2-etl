""" import_variants.py: parses SNP variants data from VCF source and loads the records to a database. """

import sys

import argparse
import logging
import sqlite3
import vcf

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
                        help='Path to the vcf file to load.')
    parser.add_argument('taxonid',
                        help='Taxon ID for the variants being loaded.')
    parser.add_argument('-c', '--create', action='store_true',
                        help='Creates new database tables (after dropping any pre-existent tables).')
    args = parser.parse_args()
    return args


def create_tables(db_conn):
    """
    DB operations to drop and then create the 'snp_variant' table.

    :param db_conn: a connection to an sqlite3 database
    :return: None
    """
    c = db_conn.cursor()
    c.execute('''DROP TABLE IF EXISTS snp_variant''')
    c.execute('''
        CREATE TABLE IF NOT EXISTS snp_variant(
          chr TEXT,
          pos INTEGER,
          id TEXT,
          ref_base TEXT,
          alt_allele TEXt,
          quality INTEGER,
          filter TEXT,
          frequency INTEGER,
          gene TEXT,
          trait_id TEXT,
          taxon_id INTEGER
        )
        ''')
    db_conn.commit()


def load_variants(db_conn, in_file, taxon_id):
    """
    Parse the VCF source and load the variants data into SNP variant table.

    :param db_conn: a connection to an sqlite3 database
    :param in_file: path to the data file
    :param taxon_id: species taxon ID for the variants being loaded
    :return: None
    """
    c = db_conn.cursor()

    var_query = """INSERT INTO snp_variant 
        (chr, pos, id, ref_base, alt_allele, quality, filter, frequency, gene, trait_id, taxon_id) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
    
    with flex_open(in_file) as in_handle:
        for rec in vcf.Reader(in_handle):
            genes = rec.INFO['CG']
            # one variant could be associated with more than one genes
            for gene in genes:
                # remove any existing 'chr' text in CHROM field
                chromosome = rec.CHROM.replace('chr', '')
                # combine alternate base(s) elements into single string
                alt_bases = ""
                if len(rec.ALT) > 0 and rec.ALT[0] is not None:
                    for base in rec.ALT:
                        alt_bases = alt_bases + base.sequence + "/"
                    alt_bases = alt_bases[:-1]

                var_record = (chromosome, rec.POS, rec.ID, rec.REF, alt_bases, rec.QUAL, rec.FILTER,
                              rec.INFO['AF'][0], gene, rec.INFO['LT'][0], taxon_id)

                try:
                    c.execute(var_query, var_record)
                except sqlite3.IntegrityError as e:
                    logger.error(f"'snp_variant' table insertion failed: {var_query}, {e}")
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

    status = load_variants(db_conn, args.filepath, args.taxonid)
    if status is True:
        print(f"  Loading variant records for <{args.taxonid}> into the database successfully completed..\n")
    else:
        print(f"  Loading variant records for <{args.taxonid}> into the database failed - check the error log..\n")


if __name__ == '__main__':
    main()
