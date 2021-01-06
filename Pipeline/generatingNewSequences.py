import argparse
import numpy as np


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
    # We initialize a dictionary, where we save for each character its number of occurrences
    counts = dict()
    # We iterate over all sequences
    for j in range(len(sequences)):
        # We access the character at position i for the current sequence j
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
        if counts[consensus] < (len(sequences) * 0.8):
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
Extract important contacts
@:param contact file with the contact information
@:param threshold threshold to characterize important contacts
@:return contacts set of residues involved in important contacts
"""


def extractImportantContacts(contact, threshold):
    contacts = set()
    i = 0

    for line in (open(contact)).readlines():
        linelist = line.split()
        if i == 0:
            i = 1
            continue
        res_1 = int(linelist[0])
        res_2 = int(linelist[1])
        conf = float(linelist[2])

        if conf < threshold:
            continue

        contacts.add(res_1)
        contacts.add(res_2)

    return contacts


"""
Extract important conserved positions.

@:param scores file with the conservation scores
@:param str_threshold threshold to characterize strictly conserved residues
@:param con_threshold threshold to characterize conserved residues
@:return set of strictly and conserved residues respectively
"""


def extractConservedPositions(scores, str_threshold, con_threshold):
    str_con = set()
    high_con = set()

    i = 0
    for line in open(scores).readlines():
        number = float(line)

        if number >= str_threshold:
            str_con.add(i)
            i += 1
            continue
        elif number >= con_threshold:
            high_con.add(i)
            i += 1
            continue
        else:
            i += 1

    return str_con, high_con


""" 
Map number to hydrophobic amino acid
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateHydrophobicAmino(number):
    switch = {0: "L", 1: "V", 2: "I", 3: "M", 4: "C", 5: "A", 6: "G", 7: "S", 8: "T", 9: "P", 10: "F", 11: "W", 12: "Y"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to polar amino acid
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translatePolarAmino(number):
    switch = {0: "E", 1: "D", 2: "N", 3: "Q", 4: "K", 5: "R", 6: "H"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to subset of amino acids
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateLeu(number):
    switch = {0: "L", 1: "V", 2: "I", 3: "M"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to subset of amino acids
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateSer(number):
    switch = {0: "S", 1: "T"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to subset of amino acids
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translatePhe(number):
    switch = {0: "F", 1: "Y", 2: "W"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to subset of amino acids
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateCha(number):
    switch = {0: "E", 1: "D", 2: "Q", 3: "N"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to subset of amino acids
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateLy(number):
    switch = {0: "K", 1: "R"}

    return switch.get(number, "Invalid Amino Acid")


""" 
Map number to amino acid
@:param number number that indicates amino acid
@:return amino acid which is mapped to number or "Invalid Amino Acid" if number is out of range
"""


def translateAmino(number):
    switch = {0: "A", 1: "R", 2: "N", 3: "D", 4: "C", 5: "E", 6: "Q", 7: "G", 8: "H", 9: "I", 10: "L", 11: "K", 12: "M",
              13: "F", 14: "P", 15: "S", 16: "T", 17: "W", 18: "Y", 19: "V"}

    return switch.get(number, "Invalid Amino Acid")


"""
Generate new protein sequence. For this a reference sequence from an alignment is used.
From this reference sequence a new protein sequence is builded according to the following scheme: 
For each residue: 
    1. If the residue is a gap in the reference sequence, it is dropped.
    2. If the residue is involved in an important contact or strictly conserved, the consensus amino acid at this 
       position is included in the new protein sequence.
    3. If the residue is conserved, the residue is substituted by a randomly chosen amino acid  which has the same 
       chemical characteristics as the original residue (can also be the original residue).
    4. If the residue is neither strictly nor conserved, a random amino acid is generated and included in the new 
       protein sequence.
@:param ref_seq reference protein sequence to take parts of natural sequence
@:param consensus consensus sequence of a alignment
@:param contacts sets of residues involved in important contacts
@:param str_con set of strictly conserved residues
@:param high_con set of conserved residues
@:return new protein sequence
"""


def generateNewSequence(ref_seq, consensus, contacts, str_con, high_con):

    """
    hydrophobic = {"L", "V", "I", "M", "C", "A", "G", "S", "T", "P", "F", "W", "Y"}
    hydrophobic_weights = [0.146, 0.10, 0.083, 0.035, 0.018, 0.136, 0.109, 0.098, 0.082, 0.072, 0.058, 0.02, 0.043]
    polar = {"E", "D", "N", "Q", "K", "R", "H"}
    polar_weights = [0.192, 0.17, 0.118, 0.117, 0.152, 0.181, 0.07]
    """

    leu = {"L", "V", "I", "M"}
    leu_weights = [0.4, 0.28, 0.22, 0.1]
    ser = {"S", "T"}
    ser_weights = [0.55, 0.45]
    phe = {"F", "Y", "W"}
    phe_weights = [0.48, 0.36, 0.16]
    cha = {"E", "D", "Q", "N"}
    cha_weights = [0.32, 0.29, 0.19, 0.2]
    ly = {"K", "R"}
    ly_weights = [0.46, 0.54]

    seq = ""

    amino_acids = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    weights = [0.092, 0.058, 0.038, 0.055, 0.013, 0.062, 0.038, 0.074, 0.022, 0.056, 0.099, 0.049, 0.027, 0.039, 0.049,
               0.067, 0.056, 0.013, 0.029, 0.064]

    for i in range(len(ref_seq)):
        if ref_seq[i] == "-" or ref_seq[i] == ".":
            continue
        else:
            if i in contacts or i in str_con:
                if consensus[i] != "-" and consensus[i] != ".":
                    seq += consensus[i]
                else:
                    seq += ref_seq[i]
            elif i in high_con:
                if ref_seq[i] in leu:
                    number = np.random.choice(len(leu), p=leu_weights)
                    amino = translateLeu(number)
                    if amino == "Invalid Amino Acid":
                        continue
                    seq += amino
                elif ref_seq[i] in ser:
                    number = np.random.choice(len(ser), p=ser_weights)
                    amino = translateSer(number)
                    if amino == "Invalid Amino Acid":
                        continue
                    seq += amino
                elif ref_seq[i] in phe:
                    number = np.random.choice(len(phe), p=phe_weights)
                    amino = translatePhe(number)
                    if amino == "Invalid Amino Acid":
                        continue
                    seq += amino
                elif ref_seq[i] in cha:
                    number = np.random.choice(len(cha), p=cha_weights)
                    amino = translateCha(number)
                    if amino == "Invalid Amino Acid":
                        continue
                    seq += amino
                elif consensus[i] in ly:
                    number = np.random.choice(len(ly), p=ly_weights)
                    amino = translateLy(number)
                    if amino == "Invalid Amino Acid":
                        continue
                    seq += amino
                else:
                    seq += ref_seq[i]
            else:
                number = np.random.choice(len(amino_acids), p=weights)
                amino = translateAmino(number)
                if amino == "Invalid Amino Acid":
                    continue
                seq += amino
            """if ref_seq[i] in hydrophobic:
                number = np.random.choice(len(hydrophobic), p=hydrophobic_weights)
                amino = translateHydrophobicAmino(number)
                seq += amino
            elif ref_seq[i] in polar:
                number = np.random.choice(len(polar), p=polar_weights)
                amino = translatePolarAmino(number)
                seq += amino
            else:
                if ref_seq[i] == " " or ref_seq[i] == "\n" or ref_seq[i] == "\t":
                    continue
                print("Undefined Amino Acid " + consensus[i])
                seq += consensus[i]
            """
    return seq


"""
Build final protein sequences. 
Enlarge the given sequence by unalignable parts from a natural sequence.
@:param indices a list of indices which map the residues of the perturbed sequence onto their position in the original
  sequence
@:param perturbed the new amino acid sequence
@:param original the original amino acid sequence
@:return the enlarged sequence
"""


def buildFinalSequence(indices, perturbed, original):
    new_seq = ""
    i = 0

    for e in range(len(indices)):
        if indices[e] != i:
            for j in range(i, indices[e]):
                new_seq += original[j]

            new_seq += perturbed[e]
            i = indices[e] + 1

    if i != len(original):
        new_seq += original[i:]
    return new_seq


"""
Maps residues from the reference sequence onto residues in the original sequence
@:param original the original sequence
@:param ref_seq the reference sequence
@:return a list of indices
"""


def getIndices(original, ref_seq):
    indices = []
    offset = 0

    while len(ref_seq) > 0:
        largest = (0, 0)
        start = 0
        number = original.find(ref_seq[0], start)
        while number != -1:
            i = 1
            length = 1
            while (number + i) < len(original) and i < len(ref_seq) and original[number + i] == ref_seq[i]:
                length += 1
                i += 1

            if length > largest[1]:
                largest = (number, length)

            start = number + length
            number = original.find(ref_seq[0], start)

        for i in range(largest[1]):
            indices.append(largest[0] + offset + i)

        ref_seq = ref_seq[0 + largest[1]:]
        original = original[largest[0] + largest[1] - 1:]
        offset += largest[0] + largest[1] - 1

    return indices


"""
Remove gaps from sequences
@:param ref_seq sequence from which gaps should be removed
@:return ref_seq without gaps
"""


def removeGaps(ref_seq):
    new_seq = ""

    for e in ref_seq:
        if e == "-" or e == ".":
            continue
        new_seq += e
    return new_seq[:-1]


"""
Find sequence to given header
@:param head header to which sequence should be found
@:param aln file of sequences to search
@:return sequence which is indicated by head or first sequence of aln if head cannot be found in aln
 (also header of this first sequence)
"""


def getSequence(head, aln):
    sequence = ""
    found = False

    for line in open(aln).readlines():
        if line.find(head) != -1:
            found = True
            continue

        if found:
            if line[0] == ">":
                break
            else:
                sequence += line[:-1]

    if not found:
        i = 0
        for line in open(aln).readlines():
            if line[0] == ">":
                if i != 0:
                    break
                else:
                    head = line
                    i = 1
            else:
                sequence += line[:-1]
    return sequence, head


def main():
    # Build an argument parser to parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=str, help='Path to the config file')
    parser.add_argument('scores', type=str, help='Path to the file with conservation scores')
    parser.add_argument('contact', type=str, help='Path to the file with contact information')
    parser.add_argument('mapping', type=str, help='Path to the output of Mapping (needed for correct order of '
                                                  'Snakemake rules)')
    parser.add_argument('output', type=str, help='Path where the output file should be written')

    args = parser.parse_args()
    scores = args.scores
    contact = args.contact
    config = args.config

    # We extract the necessary arguments
    aln = ""
    org = ""
    head = ""
    str_threshold = 0
    con_threshold = 0
    contact_threshold = 0

    for line in open(config).readlines():
        linelist = line.split("=")

        if linelist[0] == "MSA":
            aln = linelist[1][:-1]
        elif linelist[0] == "SEQUENCES":
            org = linelist[1][:-1]
        elif linelist[0] == "HEADER":
            head = linelist[1][:-1]
        elif linelist[0] == "STRICTLY_CONSERVED":
            str_threshold = float(linelist[1])
        elif linelist[0] == "CONSERVED":
            con_threshold = float(linelist[1])
        elif linelist[0] == "CORRELATED":
            contact_threshold = float(linelist[1])

    # Extract important contacts
    contacts = extractImportantContacts(contact, contact_threshold)

    # Extract conserved positions
    str_con, high_con = extractConservedPositions(scores, str_threshold, con_threshold)

    consensus = findConsensus(aln)

    # Extract reference sequence and original sequence
    # Remove gaps from reference sequence
    ref_seq, header = getSequence(head, aln)
    original, _ = getSequence(header, org)
    ref_seq = removeGaps(ref_seq)

    # Map letters of reference sequence to indices in original sequence
    indices = getIndices(original, ref_seq)

    output = args.output

    out = open(output, "w")

    for i in range(1000):
        perturbed = generateNewSequence(ref_seq, consensus, contacts, str_con, high_con)
        new_seq = buildFinalSequence(indices, perturbed, original)
        out.write(">Sequence " + str(i) + "\n")
        out.write(new_seq + "\n")

    out.close()


if __name__ == '__main__':
    main()
