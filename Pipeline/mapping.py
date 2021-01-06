import argparse
from prettytable import PrettyTable

"""
Translate one letter code to three letter code
@:param amino_acid amino acid code to translate
@:return three letter code of amino_acid or "" if not one of 20 naturally occurring amino acids
"""


def translateOneLetter(amino_acid):
    if amino_acid == "A":
        return "ALA"
    if amino_acid == "R":
        return "ARG"
    if amino_acid == "N":
        return "ASN"
    if amino_acid == "D":
        return "ASP"
    if amino_acid == "C":
        return "CYS"
    if amino_acid == "E":
        return "GLU"
    if amino_acid == "Q":
        return "GLN"
    if amino_acid == "G":
        return "GLY"
    if amino_acid == "H":
        return "HIS"
    if amino_acid == "I":
        return "ILE"
    if amino_acid == "L":
        return "LEU"
    if amino_acid == "K":
        return "LYS"
    if amino_acid == "M":
        return "MET"
    if amino_acid == "F":
        return "PHE"
    if amino_acid == "P":
        return "PRO"
    if amino_acid == "S":
        return "SER"
    if amino_acid == "T":
        return "THR"
    if amino_acid == "W":
        return "TRP"
    if amino_acid == "Y":
        return "TYR"
    if amino_acid == "V":
        return "VAL"
    return ""


"""
Merge PDB entries of two alternate residue entries
@:param first_entry first alternate residue entry
@:param second_entry second alternate residue entry
@:return merged entry which has the mean of the C-alpha coordinates of first_entry and second_entry
"""


def getNewEntry(first_entry, second_entry):
    # The amino acid of the entries without the preceding alternate indicator
    amino = first_entry[3][1:]

    # C-alpha coordinates of first entry
    first_x = float(first_entry[6])
    second_x = float(second_entry[6])
    x = (first_x + second_x) / 2

    # C-alpha coordinates of second entry
    first_y = float(first_entry[7])
    second_y = float(second_entry[7])
    y = (first_y + second_y) / 2

    first_z = float(first_entry[8])
    second_z = float(second_entry[8])

    z = (first_z + second_z) / 2

    # the newly build PDB entry with the merged coordinates
    string = first_entry[0] + " " + first_entry[1] + " " + first_entry[2] + " " + amino + " " + first_entry[4] + " "
    string += first_entry[5] + " " + str(x) + " " + str(y) + " " + str(z)
    return string


"""
Find amino acid with largest counts
@:param counts dictionary with counts for each amino acid
@:return amino acid with largest counts, gap if no amino acid has count > 0
"""


def checkPositiveAminoAcid(counts):
    maximum = -1
    new_consensus = "-"
    for key, value in counts.items():
        if key != "-" and key != ".":
            if value > maximum:
                maximum = value
                new_consensus = key
    return new_consensus


"""
Determine counts of each amino acid.
@:param i position for which counts should be determined
@:param sequences list of protein sequences
@:return dictionary of counts for amino acids
"""


def findCounts(i, sequences):
    counts = dict()

    for j in range(len(sequences)):
        c = sequences[j][i]
        if c in counts:
            counts[c] += 1
        else:
            counts[c] = 1
    return counts


"""
Determine consensus letter for a position
@:param i position for which the consensus letter should be determined
@:param sequences list of protein sequences
@:return consensus letter for position i
"""


def getConsensus(i, sequences):
    counts = findCounts(i, sequences)

    maximum = -1
    consensus = ""

    # In case of ties we take the first character found
    for key, value in counts.items():
        if value > maximum:
            maximum = value
            consensus = key

    # Check if consensus character is gap.
    # If character is gap, check if 80% or more sequence have a gap at this position.
    # If not, exchange gap by most frequent amino acid at this position.
    if consensus == "-" or consensus == ".":
        if counts[consensus] < (len(sequences) * 0.80):
            consensus = checkPositiveAminoAcid(counts)

    return consensus


"""
Determine consensus sequence
@:param aln path to the alignment
@:return the consensus sequence
"""


def findConsensus(aln):
    MSA_open = open(aln)

    sequences = []
    current_sequence = ""

    # Read sequences from alignment file
    for line in MSA_open.readlines():
        # Starting of a new header
        if line[0] == ">":
            # Append fully read sequence
            sequences.append(current_sequence)
            current_sequence = ""
        else:
            current_sequence += line[:-1]

    # Append last sequence and pop empty first sequence
    sequences.append(current_sequence)
    sequences.pop(0)

    consensus_sequence = ""

    # Determine consensus character for each position
    for i in range(len(sequences[0])):
        c = getConsensus(i, sequences)
        consensus_sequence += c

    return consensus_sequence


"""
Find sequence belonging to header
@:param header given header for which sequence should be found
@:param aln file of sequences
@:return the sequence of header or NONE if header is not found in aln
"""


def findSequence(header, aln):
    sequence = ""
    found = False

    for line in open(aln).readlines():
        if line.find(header) != -1:
            found = True
            continue

        if found:
            if line[0] == ">":
                break
            else:
                sequence += line[:-1]

    if sequence == "":
        sequence = "NONE"

    return sequence


def main():
    # We build an argument parser to parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='Path to the config file')
    parser.add_argument('conservation_scores', type=str, help='Path to the file with conservation scores')
    parser.add_argument('contact', type=str, help='Path to the file with contact information')
    parser.add_argument('output', type=str, help='Path where the output file should be written')

    arg = parser.parse_args()
    config = arg.config
    scores = arg.conservation_scores
    contact = arg.contact
    output_file = arg.output

    aln = ""
    head = ""
    pdb = ""

    for line in open(config).readlines():
        linelist = line.split("=")

        if linelist[0] == "MSA":
            aln = linelist[1][:-1]
        elif linelist[0] == "PDB":
            pdb = linelist[1][:-1]
        elif linelist[0] == "HEADER":
            head = linelist[1][:-1]

    # Find reference sequence in the file
    aligned = findSequence(head, aln)

    # If sequence cannot be found, the consensus sequence is used as alternative
    if aligned == "NONE":
        aligned = findConsensus(aln)

    # Read first L contacts (L = length of the alignment)
    contacts = dict()
    n = 0

    for line in (open(contact)).readlines():
        if n > len(aligned):
            break
        linelist = line.split()
        if len(linelist) != 2:
            continue
        res_1 = int(linelist[0])
        res_2 = int(linelist[1])

        if res_1 in contacts:
            contacts[res_1] = contacts[res_1] + 1
        else:
            contacts[res_1] = 1

        if res_2 in contacts:
            contacts[res_2] = contacts[res_2] + 1
        else:
            contacts[res_2] = 1
        n += 1

    # List of secondary structure (map secondary structure onto residue index)
    sec = [("0", "0")] * len(aligned)

    # List of PDB entries
    coordinates = []

    i = 0
    m = 1

    # In most PDB files, the first atom in the file is not the first residue of the protein sequence. To correctly
    # assign the secondary structure elements, this offset ust be tracked.
    offset = 0

    # Read secondary structures and C-alpha coordinates
    for line in open(pdb).readlines():
        linelist = line.split()
        if linelist[0] == "HELIX":
            first = int(linelist[5])
            second = int(linelist[8])
            sec[first] = (linelist[3], "HELIX")
            sec[second] = (linelist[6], "HELIX")
            for j in range(first + 1, second):
                sec[j] = ("NONE", "HELIX")
        elif linelist[0] == "SHEET":
            first = int(linelist[6])
            second = int(linelist[9])
            sec[first] = (linelist[4], "SHEET")
            sec[second] = (linelist[7], "SHEET")
            for j in range(first + 1, second):
                sec[j] = ("NONE", "SHEET")
        elif linelist[0] == "ATOM":
            if i == 0:
                offset = int(linelist[5])
                m = linelist[5]
                i = 1
            if linelist[5] != m:
                m = str(int(m) + 1)
                if linelist[2] == "CA":
                    coordinates.append(line)
            elif linelist[2] == "CA":
                coordinates.append(line)

    # Check and remove alternate PDB entries
    new_coordinates = []
    for i in range(len(coordinates)):
        line = coordinates[i]
        linelist = line.split()

        if i > 0:
            elemlist = coordinates[i - 1].split()
            if linelist[5] == elemlist[5]:
                s = getNewEntry(linelist, elemlist)
                new_coordinates.append(s)
            else:
                new_coordinates.append(line)
        else:
            new_coordinates.append(line)

    # Read conservation scores into list
    con_scores = []
    for line in open(scores).readlines():
        con_scores.append(line)

    i = 0

    # j denotes from which residue on the motif is found in the reference PDB.
    j = 0
    output = open(output_file, "w")

    table = PrettyTable()

    table.field_names = ["Alignment Position", "Position in Reference Sequence", "Amino Acid", "Secondary Structure",
                         "Coordinates", "Conservation Score", "Number of Contacts"]

    for elem in con_scores:
        if i >= len(aligned):
            break
        if aligned[i] == "-" or aligned[i] == ".":
            table.add_row([i + 1, "-", "Gap", "-", "-", elem, 0])
            i += 1
            continue
        else:
            amino = aligned[i]
            new_amino = translateOneLetter(amino)

            # Determine from which residue on the PDB sequence is represented in the alignment
            if j == 0:
                for c in range(len(new_coordinates)):
                    cor_list = (new_coordinates[c]).split()
                    if cor_list[3] == new_amino:
                        j = c
                        break

            x = new_coordinates[j]
            x_list = x.split()

            co = x_list[6] + "  " + x_list[7] + "  " + x_list[8]

            z = offset + j

            # Check for secondary structure of residue
            if sec[z] != ("0", "0"):
                planned = sec[z]
                s = planned[1]
            else:
                s = "-"

            # Check for number of contacts for residue
            try:
                con = contacts[j + 1]
            except KeyError:
                con = 0

            table.add_row([i + 1, j + 1, new_amino, s, co, elem, con])

            i += 1
            j += 1
    table_str = table.get_string()
    output.write(table_str)


if __name__ == '__main__':
    main()
