from __future__ import print_function, absolute_import
import argparse
import sqlite3

from collections import defaultdict
import re
import copy

from rdkit import Chem

from .do_fragment import parse_salt_remover
from . import command_support
from . import do_fragment
from . import fragment_algorithm
from . import environment
from . import smiles_syntax
from . import schema
from .config import DEFAULT_RULE_SELECTION_OPTIONS
from . import _compat
from . import config
import pandas as pd


class EvalError(Exception):
    pass


# =====================================
def positive_int(value):
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


# =====================================
def nonnegative_int(value):
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer or zero")
    if not (value >= 0):
        raise argparse.ArgumentTypeError("must be a positive integer or zero")
    return value


# =====================================
def open_database(dbfile):
    conn = None
    try:
        conn = sqlite3.connect(dbfile)
        return conn
    except:
        print("Problem in opening database file, exiting for now")
        exit()


# =====================================
def get_cursor(db_conn):
    return db_conn.cursor()


# =====================================
def get_fragment_options_from_args(args):
    return config.FragmentOptions(
        max_heavies=args.tmax_heavies,
        max_rotatable_bonds=args.tmax_rotatable_bonds,
        rotatable_smarts=args.trotatable_smarts,  # otherwise it's unicode
        cut_smarts=args.tcut_smarts,  # otherwise it's unicode
        num_cuts=args.tnum_cuts,
        salt_remover=args.tsalt_remover,
        method="chiral",
        min_heavies_per_const_frag=args.tmin_heavies_per_const_frag,
        min_heavies_total_const_frag=args.tmin_heavies_total_const_frag
    )


# =====================================
def _get_one_or_none(result):
    for row in result:
        return row[0]
    return None


# =====================================
def _get_all(result):
    rows = []
    for row in result:
        rows.append(row)
    return rows


# =====================================
def get_rule_smiles_id(smiles, cursor=None):
    c = cursor.execute(
        "SELECT id FROM fragment_smi WHERE fragsmi = ?", (smiles,))
    return _get_one_or_none(c)


# =====================================
def get_fingerprint_id(fingerprint, cursor=None):
    c = cursor.execute(
        "SELECT id FROM environment_fp WHERE envfp = ?", (fingerprint,))
    return _get_one_or_none(c)


# =====================================
def read_fragmentIndex_to_smiTable(conn):
    cursor = conn.cursor()
    rows = cursor.execute("SELECT id, fragsmi FROM fragment_smi")
    fragId_to_smi = {}
    for row in rows:
        fragId_to_smi[row[0]] = row[1]

    return fragId_to_smi


# =====================================
def read_envsmiId_to_envsmiTable(conn):
    cursor = conn.cursor()
    rows = cursor.execute("SELECT id, envsmi FROM environment_smi")
    envsmiId_to_smi = {}
    for row in rows:
        envsmiId_to_smi[row[0]] = row[1]

    return envsmiId_to_smi


# =====================================
def find_rule_environments_for_transform(query_pool, min_pairs=10, cursor=None, is_db_symmetric=False,
                                         db_fragId_to_fragsmi=None,
                                         db_envsmiId_to_envsmi=None):
    # make one big query can call it
    # extract (lhs_id, rhs_id, envfp_id, envsmi_id, frequency)
    sql_query = "SELECT lhs_frag, rhs_frag, envfp, envsmi, frequency, lhs_cpd, rhs_cpd, lhs_cpd_id, rhs_cpd_id FROM transformations"
    matching_rows = []
    for q in query_pool:
        subq = "(lhs_frag = %s AND envfp = %s AND frequency >= %s)" % (
            q["permuted_variable_smiles_id"], q["envfp_id"], str(min_pairs))

        if "WHERE" in sql_query:
            sql_query = sql_query + " OR " + subq
        else:
            sql_query = sql_query + " WHERE " + subq

    c = cursor.execute(sql_query)
    matching_rows = _get_all(c)

    # Now make another query and see if the query fragment exists in rhs
    # make one big query can call it
    # extract (lhs_id, rhs_id, envfp_id, envsmi_id, frequency)
    sql_query = "SELECT rhs_frag, lhs_frag, envfp, envsmi, frequency, rhs_cpd, lhs_cpd, rhs_cpd_id, lhs_cpd_id FROM transformations"
    for q in query_pool:
        subq = "(rhs_frag = %s AND envfp = %s AND frequency >= %s)" % (
            q["permuted_variable_smiles_id"], q["envfp_id"], str(min_pairs))

        if "WHERE" in sql_query:
            sql_query = sql_query + " OR " + subq
        else:
            sql_query = sql_query + " WHERE " + subq

    c = cursor.execute(sql_query)
    matching_rows = matching_rows + _get_all(c)

    # If no matches found return immediately
    if len(matching_rows) == 0:
        return []

    # Now map the lhs_id, rhs_id and envsmi_id to respective smiles
    id_to_lhsORrhs_smi = {row[0]: db_fragId_to_fragsmi[row[0]] for row in matching_rows}
    id_to_lhsORrhs_smi.update({row[1]: db_fragId_to_fragsmi[row[1]] for row in matching_rows})
    id_to_envsmi = {row[3]: db_envsmiId_to_envsmi[row[3]] for row in matching_rows}

    # group the result rows as per the query (lhs fragment id and envfp id)
    lhsId_envfpId_to_matchingRow = {}
    for row in matching_rows:

        if (row[0], row[2]) in lhsId_envfpId_to_matchingRow:
            lhsId_envfpId_to_matchingRow[(row[0], row[2])].append(row)
        else:
            lhsId_envfpId_to_matchingRow[(row[0], row[2])] = [row]

    # Final assembly
    matching_rows_output = []
    for q in query_pool:
        lhsId = q["permuted_variable_smiles_id"]
        envfpId = q["envfp_id"]

        if (lhsId, envfpId) in lhsId_envfpId_to_matchingRow:

            rows = lhsId_envfpId_to_matchingRow[(lhsId, envfpId)]
            for row in rows:
                lhs = id_to_lhsORrhs_smi[row[0]]
                rhs = id_to_lhsORrhs_smi[row[1]]
                envfp = row[2]
                envsmi = id_to_envsmi[row[3]]
                freq = row[4]
                ex_lhs_cpd_smi = row[5]
                ex_rhs_cpd_smi = row[6]
                ex_lhs_cpd_id = row[7]
                ex_rhs_cpd_id = row[8]

                matching_rows_output.append(
                    (q["frag_constant_smiles"], q["frag_variable_smiles"],
                     [lhs, rhs, envfp, envsmi, freq, ex_lhs_cpd_smi, ex_rhs_cpd_smi, ex_lhs_cpd_id, ex_rhs_cpd_id]))

    return matching_rows_output


# =====================================
class Tool(object):
    def __init__(self,
                 db_connection, fragment_options, fragment_filter):
        self.db_connection = db_connection
        self.fragment_options = fragment_options
        self.fragment_filter = fragment_filter


# =====================================
def _get_tool(klass, db_connection, args):
    fragment_options = get_fragment_options_from_args(args)
    fragment_filter = do_fragment.get_fragment_filter(fragment_options)
    return klass(
        db_connection=db_connection,
        fragment_options=fragment_options,
        fragment_filter=fragment_filter
    )


# =====================================
def get_transform_tool(db_connection, args):
    return _get_tool(TransformTool, db_connection, args)


# =====================================
class TransformTool(Tool):

    def fragment_transform_smiles(self, smiles):
        # Figure out how I'm going to fragment the input --smiles
        if "[H]" in smiles:
            # User-specified transform location
            record = do_fragment.make_hydrogen_fragment_record("query", smiles, self.fragment_filter)
        else:
            record = do_fragment.make_fragment_record_from_smiles(smiles, self.fragment_filter)
        return record

    def transform(self, transform_fragments, radius=0, min_pairs=0, min_variable_size=0,
                  max_variable_size=9999, min_constant_size=0,
                  substructure_pat=None,
                  pool=None, is_symmetric=False, db_fragId_to_fragsmi=None, db_envsmiId_to_envsmi=None):

        cursor = get_cursor(self.db_connection)
        return make_transform(
            self.db_connection, transform_fragments,
            substructure_pat=substructure_pat,
            radius=radius, min_pairs=min_pairs,
            min_variable_size=min_variable_size, max_variable_size=max_variable_size,
            min_constant_size=min_constant_size,
            pool=pool,
            cursor=cursor, is_symmetric=is_symmetric, db_fragId_to_fragsmi=db_fragId_to_fragsmi,
            db_envsmiId_to_envsmi=db_envsmiId_to_envsmi)

    def expand_variable_symmetry(self, transform_record):
        # Expand fragmentation of transform where the variable part is symmetric
        symmetry_fragments = []
        for fragment in transform_record.fragments:
            if fragment.num_cuts == 1:
                continue  # No symmetry here
            elif fragment.num_cuts == 2 and fragment.variable_symmetry_class == "11":
                if fragment.constant_symmetry_class == "11":
                    continue  # Both variable and constant are symmetric
                new_fragment = copy.copy(fragment)
                frag1, frag2 = new_fragment.constant_smiles.split(".")
                new_fragment.constant_smiles = frag2 + "." + frag1
                symmetry_fragments.append(new_fragment)

            elif fragment.num_cuts == 3 and fragment.variable_symmetry_class == '111':
                new_fragment = copy.copy(fragment)
                frag1, frag2, frag3 = new_fragment.constant_smiles.split(".")
                new_fragment.constant_smiles = frag1 + "." + frag3 + "." + frag2
                symmetry_fragments.append(new_fragment)
                new_fragment = copy.copy(fragment)
                new_fragment.constant_smiles = frag2 + "." + frag1 + "." + frag3
                symmetry_fragments.append(new_fragment)
                new_fragment = copy.copy(fragment)
                new_fragment.constant_smiles = frag2 + "." + frag3 + "." + frag1
                symmetry_fragments.append(new_fragment)
                new_fragment = copy.copy(fragment)
                new_fragment.constant_smiles = frag3 + "." + frag1 + "." + frag2
                symmetry_fragments.append(new_fragment)
                new_fragment = copy.copy(fragment)
                new_fragment.constant_smiles = frag3 + "." + frag2 + "." + frag1
                symmetry_fragments.append(new_fragment)

            elif fragment.num_cuts == 3 and fragment.variable_symmetry_class == '112':
                change_idx1, change_idx2 = int(fragment.attachment_order[0]), int(fragment.attachment_order[1])
                keep_idx = int(fragment.attachment_order[2])
                new_fragment = copy.copy(fragment)
                frags = new_fragment.constant_smiles.split(".")
                new_frags = ['', '', '']
                new_frags[keep_idx] = frags[keep_idx]
                new_frags[change_idx1] = frags[change_idx2]
                new_frags[change_idx2] = frags[change_idx1]
                new_fragment.constant_smiles = new_frags[0] + "." + new_frags[1] + "." + new_frags[2]
                symmetry_fragments.append(new_fragment)

            elif fragment.num_cuts == 3 and fragment.variable_symmetry_class == '121':
                change_idx1, change_idx2 = int(fragment.attachment_order[0]), int(fragment.attachment_order[2])
                keep_idx = int(fragment.attachment_order[1])
                new_fragment = copy.copy(fragment)
                frags = new_fragment.constant_smiles.split(".")
                new_frags = ['', '', '']
                new_frags[keep_idx] = frags[keep_idx]
                new_frags[change_idx1] = frags[change_idx2]
                new_frags[change_idx2] = frags[change_idx1]
                new_fragment.constant_smiles = new_frags[0] + "." + new_frags[1] + "." + new_frags[2]
                symmetry_fragments.append(new_fragment)

            elif fragment.num_cuts == 3 and fragment.variable_symmetry_class == '122':
                change_idx1, change_idx2 = int(fragment.attachment_order[1]), int(fragment.attachment_order[2])
                keep_idx = int(fragment.attachment_order[0])
                new_fragment = copy.copy(fragment)
                frags = new_fragment.constant_smiles.split(".")
                new_frags = ['', '', '']
                new_frags[keep_idx] = frags[keep_idx]
                new_frags[change_idx1] = frags[change_idx2]
                new_frags[change_idx2] = frags[change_idx1]
                new_fragment.constant_smiles = new_frags[0] + "." + new_frags[1] + "." + new_frags[2]
                symmetry_fragments.append(new_fragment)

        for frag in symmetry_fragments:
            transform_record.fragments.append(frag)

        return transform_record


# Enumerate all of the ways that the canonical unlabeled SMILES
# might be turned into a non-canonical labeled SMILES.

_bracket_wildcard_pat = re.compile(re.escape("[*]"))
_organic_wildcard_pat = re.compile(re.escape("*"))


def enumerate_permutations(smiles, is_symmetric=False):
    # RDKit pre-2018 used "[*]"; this changed to using a bare "*".
    if "[*]" in smiles:
        wildcard_pat = _bracket_wildcard_pat
        wildcard = "[*]"
    elif "*" in smiles:
        wildcard_pat = _organic_wildcard_pat
        wildcard = "*"

    n = smiles.count("*")
    if n == 1:
        yield "1", smiles.replace(wildcard, "[*:1]")
        return

    if n == 2:
        sub_terms = ["[*:2]", "[*:1]"]
        yield "12", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        if is_symmetric:
            return
        sub_terms = ["[*:1]", "[*:2]"]
        yield "21", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        return

    if n == 3:
        sub_terms = ["[*:3]", "[*:2]", "[*:1]"]
        yield "123", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)
        if is_symmetric:
            return

        sub_terms = ["[*:2]", "[*:3]", "[*:1]"]
        yield "132", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)

        sub_terms = ["[*:3]", "[*:1]", "[*:2]"]
        yield "213", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)

        sub_terms = ["[*:1]", "[*:3]", "[*:2]"]
        yield "231", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)

        sub_terms = ["[*:2]", "[*:1]", "[*:3]"]
        yield "312", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)

        sub_terms = ["[*:1]", "[*:2]", "[*:3]"]
        yield "321", wildcard_pat.sub(lambda pat: sub_terms.pop(), smiles)

        return

    raise AssertionError(smiles)


# The LHS only has "*", the RHS has "*:1", "*:2", ...


_weld_cache = {}


def weld_fragments(frag1, frag2):
    key = (frag1, frag2)
    value = _weld_cache.get(key, None)
    if value is not None:
        return value

    # Also cache the lhs and rhs parts because they can be reused.
    # (It's about 4% faster overall runtime on one test.)

    frag1_closures = _weld_cache.get(frag1, None)
    if frag1_closures is None:
        frag1_closures = smiles_syntax.convert_wildcards_to_closures(frag1, [1, 2, 3])
        _weld_cache[frag1] = frag1_closures

    frag2_closures = _weld_cache.get(frag2, None)
    if frag2_closures is None:
        frag2_closures = smiles_syntax.convert_labeled_wildcards_to_closures(frag2)
        _weld_cache[frag2] = frag2_closures

    welded_mol = Chem.MolFromSmiles(frag1_closures + "." + frag2_closures)
    assert welded_mol is not None, (frag1, frag2, frag1_closures + "." + frag2_closures)
    welded_smiles = Chem.MolToSmiles(welded_mol, isomericSmiles=True)

    if len(_weld_cache) > 3000:
        _weld_cache.clear()
        _weld_cache[frag1] = frag1_closures
        _weld_cache[frag2] = frag2_closures
    value = (welded_smiles, welded_mol)
    _weld_cache[key] = value
    return value


def _weld_and_filter(item):
    frag_constant_smiles, frag_variable_smiles, substructure_pat, row = item
    lhs, other_variable_smiles, envfp, envsmi, frequency, ex_lhs_cpd_smi, ex_rhs_cpd_smi, ex_lhs_cpd_id, ex_rhs_cpd_id = row
    product_smiles, new_mol = weld_fragments(frag_constant_smiles, str(other_variable_smiles))
    if substructure_pat is not None:
        # The input SMARTS can contain an explict [H],
        # which in SMARTS only matches explicit hydrogens,
        # not implicit hydrogens. It's easier to make all
        # of the hydrogens explicit than it is to adjust
        # any explicit [H] terms in the query.
        test_mol = Chem.AddHs(new_mol)
        passed_substructure_test = test_mol.HasSubstructMatch(substructure_pat)
    else:
        passed_substructure_test = True
    return (frag_constant_smiles, frag_variable_smiles, row, product_smiles, passed_substructure_test)


# ==============================================
def update_queryPool(qpool, frag_constant_smiles, frag_variable_smiles, permuted_variable_smiles,
                     permuted_variable_smiles_id, envfp_ids):
    for id in envfp_ids:
        if id is None:
            continue

        d = {"frag_constant_smiles": frag_constant_smiles, "frag_variable_smiles": frag_variable_smiles,
             "permuted_variable_smiles": permuted_variable_smiles,
             "permuted_variable_smiles_id": permuted_variable_smiles_id,
             "envfp_id": id}
        qpool.append(d)
    return qpool


# ==============================================
def make_transform(
        db_connection, transform_fragments,
        substructure_pat=None,
        radius=0, min_pairs=0, min_variable_size=0, min_constant_size=0,
        max_variable_size=9999,
        pool=None,
        cursor=None, is_symmetric=False, db_fragId_to_fragsmi=None, db_envsmiId_to_envsmi=None):
    if cursor is None:
        cursor = get_cursor(db_connection)
    assert radius in (0, 1, 2, 3, 4, 5)

    # Map from the destination SMILES to the set of rule environments
    # The RHS set contains (rule_id, rule_environment_id, is_reversed) tuples.
    output_table = {"transformed_smi": [], "constant_smi": [], "original_frag": [], "new_frag": [],
                    "envsmi": [], "rule_freq": [], "ex_lhs_cpd_smi": [], "ex_rhs_cpd_smi": [], "ex_lhs_cpd_id": [],
                    "ex_rhs_cpd_id": []}

    # Hold the welded molecules in case I need them for a substructure search

    # For each variable fragment (single, double, or triple cut) and
    # for each environment, extract all rules from the DB that start
    # with the given fragment and that has the same environment as the
    # query fragment (for all environment levels). Consider only
    # environments with radius >= the radius given as input argument.

    to_weld = []

    # A list containing all the queries. Each query is in form
    #  {frag_constant_smiles, frag_variable_smiles, permuted_variable_smiles, permuted_variable_smiles_id, envfp_id}
    query_pool = []

    # This includes the possible fragments of hydrogens
    for frag in transform_fragments:
        ## Note on terminology:
        # constant = [*]Br.[*]c1ccccc1
        # variable = c1ccc(-c2sc(-c3ccc([*])cc3)pc2[*])cc1

        # print("Processing fragment %r", frag)

        # Check if the fragmentation is allowed
        if min_variable_size and frag.variable_num_heavies < min_variable_size:
            continue

        if frag.variable_num_heavies > max_variable_size:
            continue

        if min_constant_size and frag.constant_num_heavies < min_constant_size:
            continue

        # XXX TODO: handle 'constant_with_H_smiles'?

        # In case of multiple cuts, permute the constant smiles to match the attachment order
        if frag.num_cuts > 1:
            constant_fragments = frag.constant_smiles.split(".")
            new_constant_smiles = constant_fragments[int(frag.attachment_order[0])]
            new_constant_smiles += "." + constant_fragments[int(frag.attachment_order[1])]
            if frag.num_cuts == 3:
                new_constant_smiles += "." + constant_fragments[int(frag.attachment_order[2])]
            frag.constant_smiles = new_constant_smiles

        # The variable SMILES contains unlabeled attachment points, while the
        # rule_smiles in the database contains labeled attachment points.
        # The fragment [*]CO[*] can potentially match [*:1]CO[*:2] or [*:2]CO[*:1],
        # so I need to enumerate all n! possibilities and find possible matches.

        query_possibilities = []
        for permutation, permuted_variable_smiles in enumerate_permutations(frag.variable_smiles, is_symmetric):
            permuted_variable_smiles_id = get_rule_smiles_id(permuted_variable_smiles, cursor=cursor)
            if permuted_variable_smiles_id is not None:
                query_possibilities.append((permutation, permuted_variable_smiles, permuted_variable_smiles_id))

        if not query_possibilities:
            continue

        # print(" Evaluating %d possible rule SMILES: %s",
        #      len(query_possibilities), sorted(x[0] for x in query_possibilities))

        # We now have a canonical variable part, and the assignment to the constant part.
        # Get the constant fingerprints.

        all_center_fps = environment.compute_constant_center_fingerprints_atFixRadii(
            frag.constant_smiles, radius)

        # For each possible way to represent the variable SMILES:
        #   Find all of the pairs which use the same SMILES id as the variable
        #   (The pairs are ordered so the matching SMILES is the 'from' side of the transform)
        #   The transformed SMILES goes from variable+constant -> dest_smiles+constant
        #   so weld the destination SMILES (smi2) with the constant

        for permutation, permuted_variable_smiles, permuted_variable_smiles_id in query_possibilities:
            possible_envs = environment.get_all_possible_fingerprints(
                all_center_fps, frag.variable_symmetry_class, permutation)

            envs_ids = [get_fingerprint_id(env, cursor) for env in possible_envs]

            query_pool = update_queryPool(query_pool, frag.constant_smiles, frag.variable_smiles,
                                          permuted_variable_smiles, permuted_variable_smiles_id, envs_ids)

    # Search transformation db and get rules
    # I performed bunch of queries in single call.
    # This seems to be faster, than doing each single query
    q_limit = 20
    query_pool_tmp = []
    matching_rows = []
    for control in range(len(query_pool)):
        query_pool_tmp.append(query_pool[control])

        if len(query_pool_tmp) == q_limit:
            result = find_rule_environments_for_transform(
                query_pool_tmp, min_pairs=min_pairs, cursor=cursor, db_fragId_to_fragsmi=db_fragId_to_fragsmi,
                db_envsmiId_to_envsmi=db_envsmiId_to_envsmi)
            matching_rows = matching_rows + result
            query_pool_tmp = []

    if len(query_pool_tmp) != 0:
        result = find_rule_environments_for_transform(
            query_pool_tmp, min_pairs=min_pairs, cursor=cursor, db_fragId_to_fragsmi=db_fragId_to_fragsmi,
            db_envsmiId_to_envsmi=db_envsmiId_to_envsmi)
        matching_rows = matching_rows + result

    to_weld.extend((row[0], row[1], substructure_pat, row[2]) for row in matching_rows)

    if pool is None:
        results = _compat.imap(_weld_and_filter, to_weld)
    else:
        # A chunk size of 20 seems to maximize performance.
        # Too small and there's extra pickling overhead. (Larger chunks share the same SMARTS pickle.)
        # Too large and only one process might be used for all of the welding.
        results = pool.imap(_weld_and_filter, to_weld, 20)

    for frag_constant_smiles, frag_variable_smiles, row, product_smiles, passed_substructure_test in results:
        lhs, other_variable_smiles, envfp, envsmi, frequency, ex_lhs_cpd_smi, ex_rhs_cpd_smi, ex_lhs_cpd_id, ex_rhs_cpd_id = row
        if not passed_substructure_test:
            continue

        output_table["transformed_smi"].append(product_smiles)
        output_table["constant_smi"].append(frag_constant_smiles)
        output_table["original_frag"].append(frag_variable_smiles)
        output_table["new_frag"].append(other_variable_smiles)
        output_table["envsmi"].append(envsmi)
        output_table["rule_freq"].append(frequency)
        output_table["ex_lhs_cpd_smi"].append(ex_lhs_cpd_smi)
        output_table["ex_rhs_cpd_smi"].append(ex_rhs_cpd_smi)
        output_table["ex_lhs_cpd_id"].append(ex_lhs_cpd_id)
        output_table["ex_rhs_cpd_id"].append(ex_rhs_cpd_id)

    df = pd.DataFrame.from_dict(output_table)
    return df
