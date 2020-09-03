import os
import argparse
import json
import copy

#####################################################

# This script groups the fragmentation records based on a scaffold.
# The script specifically work with the fragmentation file which is produced with mmpdb code
# The script won't work if the fragmentation file specification changes. So be careful!!

parser = argparse.ArgumentParser("grouped the fragments based on scaffold")
parser.add_argument("--fragment_folder", help="a folder containing fragmentation files generated from mmpdb algorithm")
parser.add_argument("--max_lim", help="maximum number of entries per file (if one scaffold have entries more than "
                                      "max_lim, still all entries will be write to file. "
                                      "This is only applicable when you want to put the entries from more than one "
                                      "scaffold in single file", type=int)
args = parser.parse_args()
normalizesSmiles_to_ids = {}
#####################################################
def read_frag_files(args):
    # a dictionary to hold core based fragment record
    core_to_fragmentation = {}

    # input folder
    infiles = os.listdir(args.fragment_folder)

    # Go through each file an process the data
    for file in infiles:
        infile = open(args.fragment_folder + "/" + file, "r")
        lines = infile.readlines()
        infile.close()
        for line in lines:

            if "RECORD" not in line:
                continue

            line_decoded = json.loads(line)

            # 6th entry contain's list of list. Each list is one fragmentation pattern for a given molecule
            frags = copy.deepcopy(line_decoded[5])

            # remove fragmentation from base record data
            line_decoded.pop(-1)

            # now populate the dictionary
            scaffoldIdx = 0
            for frag in frags:
                core = frag[8]
                new_record = copy.deepcopy(line_decoded)
                new_record[1] = new_record[1] + "=" + str(scaffoldIdx)
                scaffoldIdx = scaffoldIdx + 1
                new_record.append([frag])

                if core in core_to_fragmentation:
                    core_to_fragmentation[core].append(new_record)
                else:
                    core_to_fragmentation[core] = [new_record]
        print("DONE WITH FILE: ", len(core_to_fragmentation))

    print("DONE READING, TOTAL CORES: ", len(core_to_fragmentation))
    #####################################################
    # write out data
    file_counter = 0
    filew = None
    num_written = 0

    for core in core_to_fragmentation:

        # all fragmentation records for this core
        data = core_to_fragmentation[core]

        # no. of records for this core
        no_of_records = len(data)

        # if this scaffold is just appearing once, ignore it
        if no_of_records == 1:
            continue

        # If nothing is written so far. It means file is not open yet for writing. So let's open it and write out data
        if num_written == 0:
            filew = open(str(file_counter) + ".txt", "w")
            for i in range(10):
                filew.write(lines[i])
            for record in data:
                record = json.dumps(record)
                filew.write(record + "\n")

            num_written = no_of_records
            continue

        # if you are here, it means file is still open: Let's see if we can write to file or not
        if (num_written + no_of_records) > args.max_lim:

            # close the file
            filew.close()

            # increment file counter
            file_counter = file_counter + 1

            # open new file and write out data
            filew = open(str(file_counter) + ".txt", "w")
            for i in range(10):
                filew.write(lines[i])
            for record in data:
                record = json.dumps(record)
                filew.write(record + "\n")

            num_written = no_of_records
        else:
            for record in data:
                record = json.dumps(record)
                filew.write(record + "\n")
            num_written = num_written + no_of_records
    filew.close()

read_frag_files(args)
