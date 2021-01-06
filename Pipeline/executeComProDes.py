import argparse
import os
import subprocess
import time
import shutil
import multiprocessing

"""
This script builds up the necessary program calls to evaluate the newly generated sequences with ComProDes.
"""


def main():
    # Build an ArgumentParser to parse the command line arguments and add argument
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='Path to the config file')
    parser.add_argument('new_sequences', type=str, help='Path to the file with the new sequences')
    parser.add_argument('output', type=str, help='Path to file to save indices of failed sequences')

    arg = parser.parse_args()
    new_sequences = arg.new_sequences
    output = arg.output
    config = arg.config

    # Extract path to the PDB file
    pdb = ""
    for line in open(config).readlines():
        linelist = line.split("=")
        if linelist[0] == "PDB":
            pdb = linelist[1][:-1]

    # Read new sequences
    sequences = []
    current_sequence = ""

    for line in open(new_sequences).readlines():
        if line[0] == ">":
            sequences.append(current_sequence)
            current_sequence = ""
        else:
            current_sequence += line[:-1]

    # Add last sequence and remove empty first sequence.
    sequences.pop(0)
    sequences.append(current_sequence)

    # List of secondary structure (map secondary structure onto residue index)
    sec = [("0", "0")] * len(sequences[0])

    # Read secondary structures
    for line in open(pdb).readlines():
        linelist = line.split()
        if linelist[0] == "HELIX":
            first = int(linelist[5])
            second = int(linelist[8])

            # Add secondary structure to list
            # If residue is not represented by newly generated sequences, drop secondary structure
            try:
                sec[first] = (linelist[3], "H")
                sec[second] = (linelist[6], "H")
                for j in range(first + 1, second):
                    sec[j] = ("NONE", "H")
            except IndexError:
                continue
        elif linelist[0] == "SHEET":
            first = int(linelist[6])
            second = int(linelist[9])
            try:
                sec[first] = (linelist[4], "E")
                sec[second] = (linelist[7], "E")
                for j in range(first + 1, second):
                    sec[j] = ("NONE", "E")
            except IndexError:
                continue

    # String representation of secondary structures
    secondary_struc = ""

    for i in range(len(sec)):
        if sec[i] != ("0", "0"):
            planned = sec[i]
            secondary_struc += planned[1]
        else:
            secondary_struc += "C"

    # Save current working directory
    cwd = os.getcwd()

    # Prepare execution path
    path = cwd + '/ComProDes'

    # Compile script
    os.chdir(path)

    if os.system('chmod +x SEQUENCE_FITNESS.sh'):
        raise Exception('Something went wrong while compiling the evaluation script. Please check the documentation '
                        'and ensure that everything is unpacked.')

    if not os.path.exists("../../Evaluation"):
        os.mkdir("../../Evaluation")

    output_file = open("../" + output, "w")

    # Prepare number of used cores
    cores = multiprocessing.cpu_count()
    n = int(cores / 3)

    if n > 10:
        n = 10

    # Prepare sequence and evaluate
    procs = []
    for i in range(1000):
        if not os.path.exists("Sequence_" + str(i)):
            os.system("mkdir Sequence_" + str(i))
        out = open("Sequence_" + str(i) + "/Sequence_" + str(i) + ".fasta", "w")
        out.write(">Sequence_" + str(i) + ":SEQ\n")
        out.write(sequences[i] + "\n")
        out.write(">Sequence_" + str(i) + ":SEC \n")
        out.write(secondary_struc + "\n")
        out.close()
        procs.append(subprocess.Popen([path + '/SEQUENCE_FITNESS.sh', 'Sequence_' + str(i)]))

        if i % n == n - 1:
            for proc in procs:
                proc.communicate()
            procs = []

    for proc in procs:
        proc.communicate()

    for i in range(1000):
        try:
            shutil.move(path + "/Sequence_" + str(i) + "/Sequence_" + str(i) + ".SQSS", cwd[:-9] +
                        "/Evaluation/Sequence_" + str(i) + ".SQSS")
            shutil.move(path + "/Sequence_" + str(i) + "/Sequence_" + str(i) + ".DAT", cwd[:-9] +
                        "/Evaluation/Sequence_" + str(i) + ".DAT")
            os.remove(path + "/Sequence_" + str(i) + "/Sequence_" + str(i) + ".fasta")
            os.rmdir(path + "/Sequence_" + str(i))
        except IOError:
            output_file.write(str(i) + "\n")
            continue


if __name__ == '__main__':
    main()
