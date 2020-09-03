import sqlite3

# =====================================

# dbfile = "/pstore/data/cadd/awalem1/MMP_TransformationSpace/sureChEMBL/MERGINGTHEDATA/ANALYSIS/SQLITE_TABLES/transformationDBr3.sqlitdb"
# conn = sqlite3.connect(dbfile)
# cursor = conn.cursor()
# rows = cursor.execute("select * from fragment_smi")
# d = {}
# for row in rows:
#    d[row[0]] = row[1]

# print(len(d))
filer = open("/pstore/data/cadd/awalem1/MMP_TransformationSpace/sureChEMBL/MERGINGTHEDATA/ANALYSIS/testdb_toPlayAround/fragments.index", "r")
line = filer.readline()
d = {}
while line:

    if "INDEX" in line:
        line = filer.readline()
        continue

    line = line.strip()
    line = line.strip()
    line = line.split("\t")
    d[int(line[1])] = line[0]
    line = filer.readline()
print(len(d))