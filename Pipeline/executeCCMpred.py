import argparse
import os

"""
This script builds up the necessary program calls to predict the contacts of the input MSA with CCMpred.
The whole script builds up three different program calls that are executed consecutively.
The end result of the script is a list with the predicted contacts.
"""

""" 
Extract name of alignment file
@:param config path to the config file
@:return path of the MSA file
"""


def extractMSAPath(config):
    for line in open(config).readlines():
        linelist = line.split("=")
        if linelist[0] == "MSA":
            return linelist[1][:-1]


""" 
Extract name of alignment file
@:param msa_path path to the file
@:return name of the alignment file
"""


def extractMSAName(msa_path):
    msa_path_list = msa_path.split('/')

    last_elem = msa_path_list[-1]
    elem_list = last_elem.split('.')

    return elem_list[0]


""""
Returns the number of contacts that should be extracted, which is the number of possible contacts in the alignment.
@:param msa file path to the msa file
@:return number of contacts
"""

def getNumberContacts(msa):
    sequence = ""
    i = 0

    for line in open(msa).readlines():
        if line[0] == ">":
            if i == 0:
                i = 1
            else:
                break
            continue
        else:
            sequence += line[:-1]

    return len(sequence) * len(sequence)



def main():
    # Build an ArgumentParser to parse the command line arguments and add argument
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='Path to the config file')
    parser.add_argument('con_scores', type=str, help='Path to the conservation scores')
    parser.add_argument('output', type=str, help='Path to the output file')

    # Now parse the command line arguments
    arg = parser.parse_args()
    config = arg.config
    output = arg.output

    # Extract the path of the msa file
    msa = extractMSAPath(config)

    # Extract the name of the MSA file
    msa_name = extractMSAName(msa)

    # Get number of contacts that should be extracted
    n = getNumberContacts(msa)

    # First, we save the current working directory. This is important to ensure the correct location of all scripts
    # that are executed during this script.
    cwd = os.getcwd()

    # Now, we prepare the two working paths that we need for the execution.
    scripts = cwd + '/CCMpred/scripts'
    bin_path = cwd + '/CCMpred/bin'

    # In the first step of this script, the input alignment is transformed to match the CCMpred input format.
    os.chdir(scripts)

    if os.system('python3 convert_alignment.py ../../' + msa + ' fasta ../../../Results/' + msa_name + '_CCMpred_compatible.txt'):
        raise Exception('Something went wrong while preparing CCMpred input. Please check the documentation and '
                        'ensure that everything is installed correctly.')

    # Next, the main CCMpred algorithm is used to predict the contacts of the input.

    os.chdir(bin_path)
    if os.system('./ccmpred -R ../../../Results/' + msa_name + '_CCMpred_compatible.txt ../../../Results/' + msa_name
                 + '_CCMpred_contactmatrix.txt'):
        raise Exception('Something went wrong while predicting the contacts. Please check the documentation and '
                        'ensure that everything is installed correctly.')

    # Next, the contacts are extracted from the predicted contact matrix. For this first the folder needs to be
    # changed again.

    os.chdir(cwd)
    if os.system('python3 top_couplings_2.py -n ' + str(n) + ' ../Results/' + msa_name + '_CCMpred_contactmatrix.txt ' + output):
        raise Exception('Something went wrong while extracting the contacts. Please check the documentation and '
                        'ensure that everything is installed correctly.')


if __name__ == '__main__':
    main()
