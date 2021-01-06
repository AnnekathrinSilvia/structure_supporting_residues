import argparse
import os

"""
This script builds up the necessary program calls to call the C++ Evaluation script.
"""


def main():
    # Build an ArgumentParser to parse the command line arguments and add argument
    parser = argparse.ArgumentParser()
    parser.add_argument('failed', type=str, help='Path to the file with failed sequences')
    parser.add_argument('output', type=str, help='Path to the output file')

    # Now parse the command line arguments
    arg = parser.parse_args()
    output = arg.output

    # Execute the evaluation script
    os.system("./evaluatingSequences " + output)


if __name__ == '__main__':
    main()
