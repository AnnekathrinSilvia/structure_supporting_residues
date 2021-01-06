import argparse
import os
import sys

"""
This script builds up the necessary program calls to call the C++ Scoring script.
"""


def main():
    # Build an ArgumentParser to parse the command line arguments and add argument
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='Path to the config file')
    parser.add_argument('output', type=str, help='Path to the output file')

    # Now parse the command line arguments
    arg = parser.parse_args()
    config = arg.config
    output = arg.output

    # Extract the path of the MSA and the substitution matrix from the config file.
    msa = ""
    matrix = ""
    for line in open(config).readlines():
        linelist = line.split("=")
        if linelist[0] == "MSA":
            msa = linelist[1][:-1]
        elif linelist[0] == "SUBSTITUTION_MATRIX":
            matrix = linelist[1][:-1]

    if not os.path.exists("../Results"):
        os.system("mkdir ../Results")

    # Execute the scoring script
    p = os.system("./scoring " + msa + " " + matrix + " " + output)

    if p != 0:
        sys.exit("Scoring returned with exit value != 0")


if __name__ == '__main__':
    main()
