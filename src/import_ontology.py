""" import_ontology.py: parses ontology terms and mappings from OBO and GAF formatted sources and loads the records. """

import sys

import argparse
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
    # database name (required)
    parser.add_argument('database')
    parser.add_argument('-g', '--go-obo',
                        default='data/go-basic.obo.gz',
                        help="Path to the go obo file.")
    parser.add_argument('-m', '--mp-obo',
                        default='data/MPheno_OBO.ontology.gz',
                        help="Path to the mammalian phenotypes ontology file.")
    parser.add_argument('-d', '--do-obo',
                        default='data/HumanDO.obo.gz',
                        help="Path to the human disease ontology file.")
    parser.add_argument('-H', '--human-annotations',
                        default='data/Human_GeneAssociation_GO_Synteny.gaf',
                        help="Path to the file mapping human genes to GO terms.")
    parser.add_argument('-M', '--mouse-annotations',
                        default='data/Mouse_GeneAssociation_GO_Synteny.gaf',
                        help="Path to the file mapping mouse genes to GO terms.")
    parser.add_argument('-R', '--rat-annotations',
                        default='data/Rat_GeneAssociation_GO_Synteny.gaf',
                        help="Path to the file mapping mouse genes to GO terms.")
    parser.add_argument('--mouse-mp-to-gene',
                        default='data/Mouse_GeneAssociation_MP_Synteny.gaf',
                        help="Path to the file mapping mouse genes to MP terms.")
    parser.add_argument('--mouse-do-to-gene',
                        default='data/Mouse_GeneAssociation_DO_Synteny.gaf',
                        help='Path to the file mapping mouse genes to DO terms.')
    parser.add_argument('--human-do-to-gene',
                        default='data/Human_GeneAssociation_DO_Synteny.gaf',
                        help='Path to the file mapping human genes to DO terms.')
    args = parser.parse_args()
    return args


is_a = {}


def record_is_a(term, isa_term):
    """
    The GO ontology has an is_a set of relationships, which go from general
    terms to ever more specialized terms.

    To be able to handle the inheritance of terms, we need to track all.
    Later, we'll propagate all the relationships, then record all in
    the database.
    This is a bit complicated.  is_a is a transitive relationship.
        a is_a b, b is_a c => a is_a c

    Eventually, we're going to have to be able to look at c, and find all of its
    "sub-terms" (not a real term in ontology-land, I suspect).

    To make it more interesting , is_a is a many-to-many relationship:

        a is_a b
        a is_a q
        b is_a c
        d is_a c
        x is_a q

    To track this, we first build a dictionary of all the more specialized
    terms, i.e., a, b, d, and x above, with the value of those entries being
    a set of the more general terms.

    Then we're going to invert it and propagate all the relationships. In that
    way, from a parent we'll be able to find all the children, of any
    "generation".

    :param term: The term which has this is_a relationship, i.e., the more
                 specialized term.
    :param isa_term: The term in the is_a: line, the more general term.
    :return: None
    """
    if term not in is_a:
        is_a[term] = set()
    is_a[term].add(isa_term)


def save_is_a(db_con):
    c = db_con.cursor()

    # First, we need to invert the is_a dict, so that we have a dictionary of
    # general terms pointing to all their specialized terms. This is only
    # inversion.  Propagation comes later.
    inv = {}
    for specialized in is_a.keys():
        for generalized in is_a[specialized]:
            if generalized not in inv:
                inv[generalized] = set()
            inv[generalized].add(specialized)

    # Now to propagate all the specialized up to their more (and more...)
    # general terms.
    for generalized in inv.keys():
        # make a copy of the new
        s = inv[generalized]
        l = list(s)

        # Walk the list of all the more specialized terms.
        for id in l:
            # Look up the yet more specialized values, and add them to the
            # initial set (s, above), AND to the list we're walking, so we'll
            # pick up their children, too.
            try:
                # If this term has more specialized terms, process them.
                ss = inv[id]
                for s_id in ss:
                    s.add(s_id)
                    l.append(s_id)
            except:
                # If not, don't care.
                pass

    # Propagation is done. Save them in the database.
    for generalized in sorted(inv.keys()):
        for s in sorted(list(inv[generalized])):
            c.execute('''
                INSERT INTO on_pairs (
                  parent, child, relationship)
                  VALUES(?, ?, 'is_a')
            ''', (generalized, s)
            )
        # Update the generalized term in the ontology table with the length
        # of the specialized term list.
        c.execute('''
            UPDATE on_terms SET count = ? WHERE id = ?
        ''', (len(inv[generalized]), generalized))


def create_tables(db_con):
    c = db_con.cursor()

    c.execute('''DROP TABLE IF EXISTS on_terms''')
    c.execute('''
        CREATE TABLE on_terms (
          id TEXT,
          name TEXT,
          namespace TEXT,
          def TEXT,
          count INTEGER,
          PRIMARY KEY (id)
        )
    ''')
    c.execute('''CREATE INDEX on_terms_id_idx ON on_terms(id)''')
    c.execute('''CREATE INDEX on_name_idx ON on_terms(name)''')

    c.execute('''DROP TABLE IF EXISTS on_pairs''')
    c.execute('''
        CREATE TABLE on_pairs (
          parent TEXT,
          child TEXT,
          relationship TEXT
        )
    ''')
    c.execute('''CREATE INDEX rel_idx ON on_pairs(parent, relationship)''')

    c.execute('''DROP TABLE IF EXISTS gene_ontology_map''')
    c.execute('''
        CREATE TABLE gene_ontology_map (
            gene_id TEXT,
            ontology_id TEXT,
            taxonid INTEGER,
            PRIMARY KEY (gene_id, ontology_id)
        )
    ''')
    c.execute('''CREATE INDEX gene_ont_map_gene_id_idx ON
                  gene_ontology_map(gene_id, ontology_id)''')
    c.execute('''CREATE INDEX gene_ont_map_ont_id_idx ON
                  gene_ontology_map(ontology_id)''')
    c.execute('''CREATE INDEX gene_ont_map_taxonid_id_idx ON
                  gene_ontology_map(ontology_id)''')


def import_ontology(obo_file, db_con):
    c = db_con.cursor()

    first_out = False
    in_term = False
    new_term = {}
    f = flex_open(obo_file)
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line[0] == '[' and not line.startswith('[Term]'):
            in_term = False
        if line.startswith("[Term]"):
            # We're starting a new term. capture what we've seen for the current
            # one.
            # REMEMBER: We have to do this at the end of the file as well!
            if new_term:  # First time we hit [Term] we won't hve data yet.
                try:
                    c.execute('''
                        INSERT INTO on_terms (id, name, namespace, def)
                        VALUES(?, ?, ? , ?) ''',
                        (
                            new_term['id'],
                            new_term['name'],
                            new_term.get('namespace', None),
                            new_term.get('def', None),
                        )
                    )
                except sqlite3.IntegrityError:
                    print("Duplicate key!", new_term)
                    raise
            # Get ready for the next one.
            new_term = {}
            in_term = True

        # Skip if we're not in a term...
        if not in_term:
            continue

        # In a term, process it.
        if line.startswith("id: "):
            id = line.replace("id: ", "")
            new_term['id'] = id
        if line.startswith("name: "):
            name = line.replace("name: ", "")
            new_term['name'] = name
        if line.startswith("namespace: "):
            namespace = line.replace("namespace: ", "")
            new_term['namespace'] = namespace
        if line.startswith("def: "):
            definition = line.replace("def: ", "")
            new_term['def'] = definition
        if line.startswith('is_a: '):
            try:
                record_is_a(new_term['id'], line.split()[1])
            except KeyError:
                print(new_term)
                sys.exit(1)
        if line.startswith('is_obsolete: ') and \
            line.split()[1].lower() == 'true':
            # This is an obsolete term. Ignore it.
            new_term = {}    # Throw away what we've collected so far
            in_term = False  # Ignore lines until next [Term] line

    # save the last one in the DB
    c.execute('''
        INSERT INTO on_terms (id, name, namespace, def)
            VALUES(?, ?, ? , ?) ''',
                  (
                      new_term['id'],
                      new_term['name'],
                      new_term.get('namespace', None),
                      new_term.get('def', None),
                  )
              )


def import_gene_ontology_mappings(in_file, taxon_id, db_conn):
    """
    Loads the gene-ontology term mappings from GAF source in table 'gene_ontology_map'.

    :param in_file: the path to the GAF file to be loaded
    :param taxon_id: the taxon ID for the locations we are loading
    :param db_conn: a connection to an sqlite3 database
    """
    c = db_conn.cursor()

    map_query = """ INSERT OR REPLACE INTO gene_ontology_map (
                        gene_id, ontology_id, taxonid) 
                        VALUES (?, ?, ?);
                        """
    num_records = 0
    with flex_open(in_file) as in_handle:
        for line in in_handle:
            line = line.strip()
            if line[0] == "!":
                # comment line, so skip
                continue
            parts = [x.strip() for x in line.split('\t')]

            # cases when this property is: 'taxon:9606|taxon:1280'
            gene_taxonid = parts[12].replace("taxon:", "").split("|")[0]

            if int(gene_taxonid) != taxon_id:
                # species id is not matching taxon ID, so skip
                continue
            map_record = (parts[1], parts[4], gene_taxonid)
            try:
                c.execute(map_query, map_record)
                num_records = num_records + 1
            except sqlite3.IntegrityError as e:
                logger.error(f"'gene_ontology_map' table insertion Failed: {map_query}, {e}")
                db_conn.rollback()
                db_conn.close()
                return False
    return True


def main():
    args = parse_args()
    db_conn = sqlite3.connect(args.database)

    create_tables(db_conn)

    is_successful = False

    print("  Loading GO ontology terms")
    import_ontology(args.go_obo, db_conn)

    print("  Loading MP ontology terms")
    import_ontology(args.mp_obo, db_conn)

    print("  Loading DO ontology terms")
    import_ontology(args.do_obo, db_conn)

    print("  Loading GO annotations for mouse")
    is_successful = import_gene_ontology_mappings(args.mouse_annotations, 10090, db_conn)

    print("  Loading GO annotations for human")
    is_successful = import_gene_ontology_mappings(args.human_annotations, 9606, db_conn)

    print("  Loading GO annotations for rat")
    is_successful = import_gene_ontology_mappings(args.rat_annotations, 10116, db_conn)

    print("  Loading MP annotations for mouse")
    is_successful = import_gene_ontology_mappings(args.mouse_mp_to_gene, 10090, db_conn)

    print("  Loading DO annotations for mouse")
    is_successful = import_gene_ontology_mappings(args.mouse_do_to_gene, 10090, db_conn)

    print("  Loading DO annotations for human")
    is_successful = import_gene_ontology_mappings(args.human_do_to_gene, 9606, db_conn)

    save_is_a(db_conn)

    db_conn.commit()


if __name__ == '__main__':
    main()
