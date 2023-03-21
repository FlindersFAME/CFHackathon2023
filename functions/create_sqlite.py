"""
Create an SQL Lite database
"""
import gzip
import os
import sqlite3
import sys
import argparse
import time

__author__ = 'Rob Edwards'


def connect_to_db(dirname, dbname, verbose=False):
    """
    Connect to the database
    :param dirname: the database file name
    :param dbname: the name of the directory
    :param verbose: print addtional output
    :return: the database connection
    """

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
        print(f"loading roles_to_subsystems table: {dbfile} at {time.time()}", file=sys.stderr)
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
        print(f"loading mmseqs uniref table: {dbfile} at {time.time()}", file=sys.stderr)
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
        print(f"loading trmbl/sp table: {dbfile} at {time.time()}", file=sys.stderr)
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
        p = l.strip().split('\t')
        p = [x.strip() for x in p]
        acc = p[0].split('|')
        p.append(acc[1])
        p.append(acc[2])
        try:
            conn.execute("INSERT INTO trmbl_sprot VALUES (?, ?, ?, ?, ?, ?, ?, ?)", p)
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
                        "id_all": ['id', 'fig_fn', 'primary_accession', 'full_accession']}
    }

    for t in tables:
        for idx in tables[t]:
            if verbose:
                print(f"Creating index {idx} on {t} at {time.time()}", file=sys.stderr)
            conn.execute("CREATE INDEX {ix} ON {tn} ({cn})".format(ix=idx, tn=t, cn=", ".join(tables[t][idx])))
    conn.commit()

    return conn


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    connection = connect_to_db(args.o, args.d, args.v)

    subsystems(connection, '/home/edwa0468/PATRIC/subsystems-20230306.tsv.gz', args.v)
    mmseqs_uniref(connection, '/home/edwa0468/UniRef/mmseqs_uniref.tsv.gz', args.v)
    trembl_sprot(connection, '/home/edwa0468/UniRef/trembl.proc.gz', args.v)
    trembl_sprot(connection, '/home/edwa0468/UniRef/trembl.proc.gz', args.v)
    create_indices(connection, args.v)
