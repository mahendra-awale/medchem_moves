from rdkit import Chem
import subprocess
import sqlite3

code_base = "C:\\Users\\awalem1\\Desktop\\MMPDB_cpdGenerator_ForWebS\mmpdb"

# directory where transformation database are stored
db = 'C:\\Users\\awalem1\\Documents\\rocheDB_r3.tdb'

##############################################
def decode_output(output_string):
    output_string = output_string.decode()
    output_string = output_string.replace("*:1", "R1")
    output_string = output_string.replace("*:2", "R1")
    output_string = output_string.replace("*:3", "R1")
    output_string = output_string.replace("*", "[R1]")

    if output_string == "None":
        return None
    else:
        output = []
        output_lines = output_string.split("\n")
        filew = open("C:\\Users\\awalem1\\Desktop\\tmp.txt", "w")
        for line in output_lines:
            filew.write(line+"\n")
        filew.close()
        print(output_lines)
##############################################
def generate_cpds():

    # run the compound generator
    p = subprocess.run(
        ["python", code_base,
         "mmpCompoundGenerator",
         "--transformation_db", db,
         "--tmin-pairs", str(1),
         '--tsmiles', "Clc1ccc(cc1)C(c2ccccc2)N3CCN(CC3)CCOCC(=O)O",
         '--replaceGroup', "None",
         "--tradius", str(3),
         "--tjobs", "8"],
        stdout=subprocess.PIPE, shell=False)

    # output (its in byte format)
    output = p.stdout
    output = decode_output(output)
    return output

generate_cpds()