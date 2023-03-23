"""
Convert functions to subsystems
"""

import os
import sys
import argparse
import sqlite3

__author__ = 'Rob Edwards'


def connect_to_db(db, verbose=False):
    """
    Connect to the database
    :param db: the path and name of the database
    :param verbose: print addtional output
    :return: the database connection
    """

    if not os.path.exists(db):
        print(f"Database {db} does not exist", file=sys.stderr)
        sys.exit(0)

    try:
        connx = sqlite3.connect(db)
    except sqlite3.Error as e:
        sys.stderr.write(f"ERROR connecting to database: {db}\n")
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


def search_file(conn, file, column, verbose=False):
    with open(file, 'r') as f:
        for l in f:
            query = l.strip()
            if "\t" in l:
                p = l.split("\t")
                query = p[column]
            if verbose:
                print(l, file=sys.stderr)
            cur = conn.cursor()
            try:
                conn.execute("select * from roles_to_subsystems where role = ?", [query])
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write(f"\nWhile insert on: {p}\n")
                sys.exit()
            p = cur.fetchone()

            if p:
                print("\t".join(p))


if __name__ == "__main__":
    defaultdb = '/home/edwa0468/CF_Hackathon2023/database/functions.sqlite'
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--file', help='file of seed/fig/patric functions, one per line', required=True)
    parser.add_argument('-c', '--column', help='column from tab-separated (default=0)', default=0, type=int)
    parser.add_argument('-d', '--db', help=f'database (incl. path) {defaultdb}', default=defaultdb)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    conx = connect_to_db(args.db, args.verbose)
    search_file(conx, args.file, args.column, args.verbose)
    disconnect(conx, args.verbose)
