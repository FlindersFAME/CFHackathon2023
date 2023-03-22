"""
Create an SQL Lite database
"""
import gzip
import os
import sqlite3
import sys
import argparse
from datetime import datetime

__author__ = 'Rob Edwards'


def connect_to_db(dbname, dirname, verbose=False):
    """
    Connect to the database
    :param dbname: the name of the directory
    :param dirname: the database file name
    :param verbose: print addtional output
    :return: the database connection
    """


    os.makedirs(dirname, exist_ok=True)

    try:
        if verbose:
            sys.stderr.write(f"Connecting to {os.path.join(dirname, dbname)}\n")
        connx = sqlite3.connect(os.path.join(dirname, dbname))
    except sqlite3.Error as e:
        sys.stderr.write(f"ERROR Creating database: {os.path.join(dirname, dbname)}\n")
        sys.stderr.write(e)
        sys.exit(-1)

    if verbose:
        sys.stderr.write(f"Connected to database: {sqlite3.version}\n")

    return connx


def disconnect(conn, verbose=False):
    """
    Disconnect the database and ensure we've saved all changes
    :param conn: the database connection
    :param verbose: print addtional output
    :return:
    """

    if conn:
        conn.commit()
        conn.close()
    elif verbose:
        sys.stderr.write("There was no database connection!\n")


def subsystems(conn, dbfile='', verbose=False):
    """
    create a subsystems database

    :param conn: the database connection
    :param dbfile: the database file name
    :param verbose: more output
    """

    if verbose:
        print(f"loading roles_to_subsystems table: {dbfile} at {datetime.now()}", file=sys.stderr)
    if not os.path.exists(dbfile):
        sys.stderr.write(f"ERROR: {dbfile} does not exist\n")
        sys.exit(-1)
    conn.execute("""
        CREATE TABLE roles_to_subsystems (
        role_id INTEGER PRIMARY KEY, 
        role TEXT, 
        subsystem TEXT,
        subclass TEXT,
        class TEXT,
        superclass TEXT
        )
        """)
    conn.commit()

    if dbfile.endswith('.gz'):
        fin = gzip.open(dbfile, 'rt')
    else:
        fin = open(dbfile, 'r')

    for l in fin:
        p = l.strip().split('\t')
        p = [x.strip() for x in p]
        try:
            conn.execute("INSERT INTO roles_to_subsystems (role, subsystem, subclass, class, superclass) "
                         "VALUES (?, ?, ?, ?, ?)", p)
        except sqlite3.OperationalError as e:
            sys.stderr.write("{}".format(e))
            sys.stderr.write(f"\nWhile insert on: {p}\n")
            sys.exit()
    conn.commit()


def mmseqs_uniref(conn, dbfile, verbose=False):
    """
    create a subsystems database

    :param conn: the database connection
    :param dbfile: the databae file
    :param verbose: more output
    """

    if verbose:
        print(f"loading mmseqs uniref table: {dbfile} at {datetime.now()}", file=sys.stderr)
    if not os.path.exists(dbfile):
        sys.stderr.write(f"ERROR: {dbfile} does not exist\n")
        sys.exit(-1)
    conn.execute("""
        CREATE TABLE mmseqs_uniref (
        uniref_id TEXT PRIMARY KEY, 
        uniref_fn TEXT,
        num_proteins INTEGER,
        taxonomy TEXT,
        taxid INTEGER,
        repID TEXT
        )
        """)
    conn.commit()

    if dbfile.endswith('.gz'):
        fin = gzip.open(dbfile, 'rt')
    else:
        fin = open(dbfile, 'r')

    for l in fin:
        p = l.strip().split('\t')
        p = [x.strip() for x in p]
        p[2] = int(p[2].replace('n=', ''))
        p[3] = p[3].replace('Tax=', '')
        p[4] = int(p[4].replace('TaxID=', ''))
        p[5] = p[5].replace('RepID=', '')
        try:
            conn.execute("INSERT INTO mmseqs_uniref VALUES (?, ?, ?, ?, ?, ?)", p)
        except sqlite3.OperationalError as e:
            sys.stderr.write("{}".format(e))
            sys.stderr.write(f"\nWhile insert on: {p}\n")
            sys.exit()
    conn.commit()


def trembl_sprot(conn, dbfile, verbose=False):
    """
    Create a table for the trmbl data. This has additional identifier columns
    """

    if verbose:
        print(f"loading trmbl/sp table: {dbfile} at {datetime.now()}", file=sys.stderr)
    if not os.path.exists(dbfile):
        sys.stderr.write(f"ERROR: {dbfile} does not exist\n")
        sys.exit(-1)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS trmbl_sprot (
        id TEXT PRIMARY KEY, 
        md5 TEXT,
        pgfam TEXT,
        fig_fn TEXT,
        uniref_fn TEXT,
        kv_pairs TEXT,
        primary_accession TEXT,
        full_accession TEXT
        )
        """)
    conn.commit()

    if dbfile.endswith('.gz'):
        fin = gzip.open(dbfile, 'rt')
    else:
        fin = open(dbfile, 'r')

    for l in fin:
        if l.startswith('id'):
            continue
        p = l.strip().split('\t')
        p = [x.strip() for x in p]
        acc = p[0].split('|')
        if len(acc) > 2:
            p.append(acc[1])
            p.append(acc[2])
        else:
            p.append('')
            p.append('')
        try:
            conn.execute("INSERT INTO trmbl_sprot VALUES (?, ?, ?, ?, ?, ?, ?, ?)", p)
        except sqlite3.OperationalError as e:
            sys.stderr.write("{}".format(e))
            sys.stderr.write(f"\nWhile insert on: {p}\n")
            sys.exit()
    conn.commit()


def id_mapping(conn, dbfile, verbose=False):
    """
    Create a table for the trmbl data. This has additional identifier columns
    """

    if verbose:
        print(f"loading idmapping table: {dbfile} at {datetime.now()}", file=sys.stderr)
    if not os.path.exists(dbfile):
        sys.stderr.write(f"ERROR: {dbfile} does not exist\n")
        sys.exit(-1)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS id_map (
            UniProtKB_AC TEXT,
            UniProtKB_ID TEXT,
            GeneID TEXT,
            RefSeq TEXT,
            GI TEXT,
            PDB TEXT,
            GO TEXT,
            UniRef100 TEXT,
            UniRef90 TEXT,
            UniRef50 TEXT,
            UniParc TEXT,
            PIR TEXT,
            NCBI_taxon TEXT,
            MIM TEXT,
            UniGene TEXT,
            PubMed TEXT,
            EMBL TEXT,
            EMBL_CDS TEXT,
            Ensembl TEXT,
            Ensembl_TRS TEXT,
            Ensembl_PRO TEXT,
            Additional_PubMe TEXT,
        )
        """)
    conn.commit()

    if dbfile.endswith('.gz'):
        fin = gzip.open(dbfile, 'rt')
    else:
        fin = open(dbfile, 'r')

    for l in fin:
        if l.startswith('id'):
            continue
        p = l.strip().split('\t')
        p = [x.strip() for x in p]

        try:
            conn.execute("INSERT INTO trmbl_sprot VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", p)
        except sqlite3.OperationalError as e:
            sys.stderr.write("{}".format(e))
            sys.stderr.write(f"\nWhile insert on: {p}\n")
            sys.exit()
    conn.commit()




def create_indices(conn, verbose=False):
    """
    Create some useful indices. Note that the PRIMARY KEY columns are indexed by default!
    :param conn: The database connection
    :param verbose: print addtional output
    :return:
    """

    if verbose:
        sys.stderr.write("Creating indices\n")

    tables = {
        "roles_to_subsystems": {"ss_role": ["role", "subsystem"],
                                "ss_role_class": ["role", "subsystem", "subclass", "class", "superclass"]},
        "mmseqs_uniref": {
            "uniref": ["uniref_id", "repID"],
            "unireffn": ["uniref_id", "uniref_fn"]
        },
        "trmbl_sprot": {"id_fig": ['id', 'fig_fn'],
                        "id_pa": ['id', 'primary_accession'],
                        "id_fa": ['id', 'full_accession'],
                        "id_all": ['id', 'fig_fn', 'primary_accession', 'full_accession']},
        "id_map": {
            "uniuni": ['uniparc', 'uniref50']
        }
    }

    for t in tables:
        for idx in tables[t]:
            if verbose:
                print(f"Creating index {idx} on {t} at {datetime.now()}", file=sys.stderr)
            conn.execute("CREATE INDEX IF NOT EXISTS {ix} ON {tn} ({cn})".format(ix=idx, tn=t, cn=", ".join(tables[t][idx])))
    conn.commit()

    return conn


if __name__ == "__main__":

    databases = {
        'subsystems' : '/home/edwa0468/PATRIC/subsystems-20230306.tsv.gz',
        'mmseqs_uniref' : '/home/edwa0468/UniRef/mmseqs_uniref.tsv.gz',
        'sprot': '/home/edwa0468/UniRef/sprot.proc.gz',
        'trembl': '/home/edwa0468/UniRef/trembl.proc.gz',
        'idmapping': '/home/edwa0468/UniRef/idmapping_selected.tab.gz'
    }

    parser = argparse.ArgumentParser(description='Create an SQLlite database and load databases')
    parser.add_argument('-f', '--dbfile', help=f'database file', required=True)
    parser.add_argument('-d', '--directory', help='database directory', required=True)
    grp = parser.add_argument_group('database to load')
    grp.add_argument('-a', '--all', help='load all databases', action='store_true')
    grp.add_argument('-s', '--subsystems', help=f"Subsystems: {databases['subsystems']}", action='store_true')
    grp.add_argument('-m', '--mmseqs', help=f"mmseqs_uniref: {databases['mmseqs_uniref']}", action='store_true')
    grp.add_argument('-r', '--sprot', help=f"sprot: {databases['sprot']}", action='store_true')
    grp.add_argument('-t', '--trembl', help=f"trembl: {databases['trembl']}", action='store_true')
    grp.add_argument('-i', '--ids', help=f"Subsystems: {databases['idmapping']}", action='store_true')

    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    connection = connect_to_db(args.dbfile, args.directory, args.verbose)

    if args.all or args.subsystems:
        subsystems(connection, databases['subsystems'], args.verbose)
    if args.all or args.mmseqs:
        mmseqs_uniref(connection, databases['mmseqs_uniref'], args.verbose)
    if args.all or args.sprot:
        trembl_sprot(connection, databases['sprot'], args.verbose)
    if args.all or args.trembl:
        trembl_sprot(connection, databases['trembl'], args.verbose)
    if args.all or args.ids:
        id_mapping(connection, databases['idmapping'], args.verbose)

    # create indices as this only works with new indices
    create_indices(connection, args.verbose)
