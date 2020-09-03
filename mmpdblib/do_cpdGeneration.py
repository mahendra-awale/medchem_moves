# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from __future__ import print_function

import sys
import time
import multiprocessing
import numpy as np
from rdkit import Chem

from . import command_support
from . import dbutils
# from . import cpdGeneration_algorithm
from . import cpdGeneration_algorithm_optimized

from . import fileio
# from .cpdGeneration_algorithm import open_database
from .cpdGeneration_algorithm_optimized import open_database

import pandas as pd
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

########################
def add_diff_heavies(dataFrame):
    '''
    Compute the difference between new and old fragment
    and update the dataframe
    '''
    diff = []
    for idx, row in dataFrame.iterrows():
        f1 = row.original_frag
        f2 = row.new_frag
        f1 = Chem.MolFromSmiles(f1)
        f2 = Chem.MolFromSmiles(f2)

        hac1 = 0
        hac2 = 0
        if f1 is not None:
            hac1 = f1.GetNumHeavyAtoms()
        if f2 is not None:
            hac1 = f2.GetNumHeavyAtoms()
        diff.append(np.abs(hac1 - hac2))

    dataFrame.insert(len(dataFrame.columns), "heavies_diff", diff)
    return dataFrame

def remove_star_atom(mol):
    for atom in mol.GetAtoms():
         if atom.GetSymbol() == "*" or atom.GetSymbol() == "[*]":
             atom.SetAtomicNum(2)
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[He]'))
    return mol


def is_containSubstructure(mol, patt):
    if mol is None or patt is None:
        return True
    mol = remove_star_atom(mol)
    patt = remove_star_atom(patt)
    return mol.HasSubstructMatch(patt)
    

########################
def get_fragmentAtomMapping(replaceGroup, qmolecule):
    '''
    Mapped the replaceGroup to query molecule and returns
    corresponding atom indexes from query molecule
    replaceGroup: smiles notation of replaceGroup
    qmolecule: smiles notation of complete query molecule
    '''

    replaceMol = Chem.MolFromSmarts(replaceGroup)
    qmol = Chem.MolFromSmiles(qmolecule, sanitize=False)
    matches = list(qmol.GetSubstructMatches(replaceMol))

    if len(matches) == 0:
        return matches

    # ommit the * atom index from list
    # logic is that: atom mapping order result from qmol
    # equals to order of atoms in replaceGroup. Which means:
    # if * atom is having index 5 in replaceGroup
    # the corresponding matching atom from query molecule will be
    # at 6th position in result list.
    filter_matches = []
    for mat in matches:
        for i in range(len(mat)):
            if replaceMol.GetAtomWithIdx(i).GetSymbol() == "*":
                continue
            else:
                filter_matches.append(mat[i])
    return filter_matches


########################
# Helper function to make a new function which format time deltas
# so all the "."s are lined up and they use the minimal amount of
# left-padding.
def get_time_delta_formatter(max_dt):
    s = "%.1f" % (max_dt,)
    num_digits = len(s)
    fmt = "%" + str(num_digits) + ".1f"

    def format_dt(dt):
        return fmt % (dt,)

    return format_dt


########################
def mmpCompoundGenerator(parser, args):
    # Generator
    radius = args.tradius
    assert radius in list("012345"), radius
    radius = int(radius)
    min_pairs = int(args.tmin_pairs)
    min_variable_size = args.tmin_variable_size
    max_variable_size = args.tmax_variable_size
    assert max_variable_size > min_variable_size, "max-variable-size must be greater than min-variable-size"
    min_constant_size = args.tmin_constant_size

    start_time = time.time()
    db_connection = open_database(args.transformation_db)
    open_time = time.time()

    # I preferred to do this, as querying several times this database is seems to be slow
    # If you make one big query and execute it, then you will encounter other sqlite runtime error
    # for instance max. tree dept limit
    # read fragment index to fragment smi table
    # Off course there is work around, but its like more coding for low performance gain
    fragId_to_fragsmi = cpdGeneration_algorithm_optimized.read_fragmentIndex_to_smiTable(db_connection)

    # read envsmi index to envsmi table
    envsmiId_to_envsmi = cpdGeneration_algorithm_optimized.read_envsmiId_to_envsmiTable(db_connection)

    if args.tsubstructure:
        substructure_pat = Chem.MolFromSmarts(args.tsubstructure)
        if substructure_pat is None:
            parser.error("Cannot parse --substructure %r" % (args.tsubstructure,))
    else:
        substructure_pat = None

    # get the transform tool
    transform_tool = cpdGeneration_algorithm_optimized.get_transform_tool(db_connection, args)

    # output data-frame
    output_df = None

    # query smiles
    smi = args.tsmiles

    # fragment entire compound
    transform_record = transform_tool.fragment_transform_smiles(smi)
    transform_record = transform_tool.expand_variable_symmetry(transform_record)

    if transform_record.errmsg:
        sys.stdout.write("ERROR: Unable to fragment --smiles %r: %s" % (args.tsmiles, transform_record.errmsg))
        exit()

    # Check if the replace-group fragment specified by users is
    # present in fragmentation of entire molecule or not
    replaceGroup_Found = False
    replaceGroup_Mol = None
    fragments = None

    # Three possible scenes for fragmentation's:
    # 1) if replace group is not provided i.e. None then consider all fragmentation's
    # 2) if replace group is single atom then consider that specific fragmentation
    # 3) if replace group is more than one atom, then re-fragment replaceGroup and
    # consider all those fragments as queries
    possible_frags = {}
    original_constatPart_as_mol = None
    if args.tconstant_smi is not None:
        original_constatPart_as_mol = Chem.MolFromSmiles(args.tconstant_smi, sanitize=True)
 
    if args.replaceGroup is not None and args.replaceGroup != "None":
        replaceGroup_Mol = Chem.MolFromSmiles(args.replaceGroup, sanitize=True)
        args.replaceGroup = Chem.MolToSmiles(replaceGroup_Mol, isomericSmiles=True)
        fragments = []
        # * and [*] both are same, but depending upon rdkit version you might get * or [*]
        # For comparision here just make everything *
        replaceGroup = args.replaceGroup.replace("[*]", "*")
        for fragment_record in transform_record.fragments:
            variable_smiles = fragment_record.variable_smiles.replace("[*]", "*")
            if replaceGroup == variable_smiles:
                replaceGroup_Found = True
                possible_frags[replaceGroup] = replaceGroup
                fragments.append(fragment_record)

        # if replaceGroup not found, then exit right away
        if not replaceGroup_Found:
            sys.stdout.write("ERROR: Replace group %s not found in fragmentation of input molecule %s" % (
                replaceGroup, args.tsmiles))
            exit()

        # if replace group is not a single atom group, then do fragmentation of replaceGroup
        if replaceGroup not in ["*Cl", "*F", "I", "*Br", "*C", "*c", "*N", "*n", "*O", "*o", "*S", "*s", "*P", "*[H]"]:
   
            # get replaceGroup atom mapping in query molecule
            #replaceGroup_atomMapping = get_fragmentAtomMapping(replaceGroup, args.tsmiles)

            # replace * in replaceGroup by Argon (its dummy atom)
            for i in range(len(replaceGroup_Mol.GetAtoms())):
                if replaceGroup_Mol.GetAtomWithIdx(i).GetSymbol() == "*":
                    replaceGroup_Mol.GetAtomWithIdx(i).SetAtomicNum(18)

            # do the fragmentation of replaceGroup molecule
            replaceGroup_smi = Chem.MolToSmiles(replaceGroup_Mol, isomericSmiles=True)
            replaceGrp_transform_record = transform_tool.fragment_transform_smiles(replaceGroup_smi)
            replaceGrp_transform_record = transform_tool.expand_variable_symmetry(replaceGrp_transform_record)
            for fragment_record in replaceGrp_transform_record.fragments:
                variable_smiles = fragment_record.variable_smiles.replace("[*]", "*")
                variable_smiles = variable_smiles.replace("[Ar]", "*")
                possible_frags[variable_smiles] = variable_smiles
            
            fragments = []
            for frag_record in transform_record.fragments:
                variable_smiles = frag_record.variable_smiles.replace("[*]", "*")
                if variable_smiles in possible_frags:
                    fragments.append(frag_record) 
    else:
        fragments = transform_record.fragments

    # filter the fragments to restrict to the ones, whose constant part is present in original constant part
    fragments_filter = []
    for frag in fragments:
        if original_constatPart_as_mol is not None and args.tconstant_smi.count("*") == 1:
            constantPart = Chem.MolFromSmiles(frag.constant_smiles, sanitize=True)
            is_pass = is_containSubstructure(constantPart, original_constatPart_as_mol)
            if is_pass:
                fragments_filter.append(frag)
            else:
                continue
        else:
            fragments_filter.append(frag)


    try:
        pool = multiprocessing.Pool(processes=args.tjobs)
        output_df = transform_tool.transform(fragments_filter,
                                             radius=radius,
                                             min_pairs=min_pairs,
                                             min_variable_size=min_variable_size,
                                             max_variable_size=max_variable_size,
                                             min_constant_size=min_constant_size,
                                             substructure_pat=substructure_pat,
                                             pool=pool, db_fragId_to_fragsmi=fragId_to_fragsmi,
                                             db_envsmiId_to_envsmi=envsmiId_to_envsmi
                                             )

        # add the original smiles to data-frame
        smis = [smi] * len(output_df)
        output_df.insert(0, "original_smi", smis)
        
        if len(output_df) == 0:
            sys.stdout.write("ERROR: Everything was good, but no rules were found")
            exit()

        # sort based on freq
        output_df = output_df.sort_values(["rule_freq"], ascending=False)
        output_df.reset_index(inplace=True, drop=True)

    except cpdGeneration_algorithm_optimized.EvalError as err:
        sys.stdout.write("ERROR: %s\nExiting.\n" % (err,))
        exit()

    # add diff heavies column
    output_df = add_diff_heavies(output_df)
    
    # remove duplicates
    output_df = output_df.drop_duplicates(["transformed_smi"])
    output_df.reset_index(inplace=True, drop=True)

    if args.toutput is not None:
        output_df.to_csv(args.toutput, sep="\t", index=False)
        sys.stdout.write("DONE: output written to %s file" %args.toutput)
        exit()

    columns = output_df.columns
    output_str = "\t".join(output_df.columns)
    for idx, row in output_df.iterrows():
        output_str = output_str + "\n" + "\t".join([str(row[col]) for col in columns])
    sys.stdout.write(output_str)

    if args.toutput is not None:
        output_df.to_csv(args.toutput, sep="\t", index=False)
########################
