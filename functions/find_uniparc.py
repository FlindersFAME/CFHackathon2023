"""
Search for uniparc IDs in the sqlite database
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


def search_file(conn, file, verbose=False):
    with open(file, 'r') as f:
        for l in f:
            l = l.strip()
            l = l.replace('UniRef50_', '')
            cur = conn.cursor()
            try:
                conn.execute("select * from id_map where uniparc = ?", l)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write(f"\nWhile insert on: {p}\n")
                sys.exit()
            p = cur.fetchone()

            if p:
                t = "\t".join(p)
                print(f"{l}\t{t}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--file', help='file of UniParc IDs to test', required=True)
    parser.add_argument('-d', '--db', help='database (incl. path)', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    conx = connect_to_db(args.db, args.verbose)
    search_file(conx, args.file, args.verbose)
    disconnect(conx, args.verbose)
