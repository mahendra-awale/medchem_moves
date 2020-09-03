import sqlite3
import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--mmpdb", help="a mmpdb sql database as input", type=str)
parser.add_argument("--output_file", help="a name of output file", type=str)
args = parser.parse_args()


#########################################################
def open_sqlitedb_reader(mmpdb_filename):
    """
    opens up the connection to matched molecular series database
    :param mmsdb_filename: path to mmsdb
    :return: connection to matched molecular series database
    """
    if not os.path.isfile(mmpdb_filename):
        print("error, database file %s does not exists" % mmpdb_filename)
        return None

    db = sqlite3.connect(mmpdb_filename)
    conn = db.cursor()
    conn.execute("PRAGMA case_sensitive_like = true")
    return conn


#########################################################
def get_publicCpdIdTable(conn):
    privateCpdid_to_publicCpdid = {}

    sql = "SELECT * from compound"
    result = conn.execute(sql)
    row_header_description = result.description
    row_header = [x[0] for x in row_header_description]
    rows = result.fetchall()

    for r in rows:
        data_dict = dict(zip(row_header, r))
        priv_id = data_dict["id"]
        pub_id = data_dict["public_id"]
        privateCpdid_to_publicCpdid[priv_id] = pub_id

    return privateCpdid_to_publicCpdid

#########################################################
def get_idToRule_Smiles(conn):
    id_to_rule_smiles = {}

    sql = "SELECT * from rule_smiles"
    result = conn.execute(sql)
    row_header_description = result.description
    row_header = [x[0] for x in row_header_description]
    rows = result.fetchall()

    for r in rows:
        data_dict = dict(zip(row_header, r))
        id = data_dict["id"]
        smiles = data_dict["smiles"]
        id_to_rule_smiles[id] = smiles

    return id_to_rule_smiles


#########################################################
def process(args):
    # open mmpdb reader
    conn = open_sqlitedb_reader(args.mmpdb)

    # Get private compound id to public compound id map
    privateCpdid_to_publicCpdid = get_publicCpdIdTable(conn)

    # Get id to variable smiles map
    id_to_ruleSmiles = get_idToRule_Smiles(conn)

    # convert it to flight file
    sql = "SELECT * from pair " \
          "INNER JOIN rule_environment ON rule_environment.id = pair.rule_environment_id " \
          "INNER JOIN rule ON rule_environment.rule_id = rule.id " \
          "INNER JOIN environment_fingerprint ON environment_fingerprint.id = rule_environment.environment_fingerprint_id " \
          "INNER JOIN constant_smiles ON pair.constant_id = constant_smiles.id"

    result = conn.execute(sql)
    row_header_description = result.description
    row_header = [x[0] for x in row_header_description]
    rows = result.fetchall()

    output = {"from_cpdid": [], "tocpd_id": [], "constant": [], "rule": [], "radius": [], "fp": []}
    for r in rows:
        data_dict = dict(zip(row_header, r))
        cpd1 = privateCpdid_to_publicCpdid[data_dict["compound1_id"]]
        cpd2 = privateCpdid_to_publicCpdid[data_dict["compound2_id"]]

        # check for compound identity (= is added to compound id in previous step while indexing for scaffold)
        if cpd1.split("=")[0] == cpd2.split("=")[0]:
            continue

        constant = data_dict["smiles"]
        radius = data_dict["radius"]
        fp = data_dict["fingerprint"]
        from_smiles = id_to_ruleSmiles[data_dict["from_smiles_id"]]
        to_smiles = id_to_ruleSmiles[data_dict["to_smiles_id"]]

        output["from_cpdid"].append(cpd1)
        output["tocpd_id"].append(cpd2)
        output["constant"].append(constant)
        output["rule"].append(from_smiles + ">>" + to_smiles)
        output["radius"].append(radius)
        output["fp"].append(fp)

    cols = list(output.keys())
    df = pd.DataFrame.from_dict(output)
    df.to_csv(args.output_file, sep="\t", columns=cols, index=False)


# process the database
process(args)
