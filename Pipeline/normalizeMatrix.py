import argparse
import os

"""
This script normalizes a given substitution matrix. 
To normalize the matrix the diagonal values are substituted by the mean value of the diagonal values.
"""

"""
Calculate the mean value of the diagonal entries.
@:param matrix given matrix
@:param isInt boolean value indicating if matrix contains only integer values
@:return mean value of diagonal values
"""


def getMean(matrix, isInt):
    sum = 0

    for n in range(len(matrix)):
        for m in range(len(matrix[n])):
            if n == m:
                sum += matrix[n][m]

    if isInt:
        sum = int(sum / len(matrix))
    else:
        sum = sum / len(matrix)

    return sum


def main():
    # Build an ArgumentParser to parse the command line arguments and add argument
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix', type=str, help='Path to the unnormalized Matrix file')

    # Now parse the command line arguments
    arg = parser.parse_args()
    mat = arg.matrix

    output = ""

    matrixlist = mat.split(".")
    output = matrixlist[0] + "_normalized.txt"

    # Now we parse the substitution matrix into a list of list.
    i = 0
    matrix = []
    amino_acids = []
    isInt = False

    for line in (open(mat)).readlines():
        if i == 0:
            i = 1
            continue
        else:
            linelist = line.split()
            amino = linelist[0]
            amino_acids.append(amino)
            linelist = linelist[1:]
            current_list = []

            for n in range(len(linelist)):
                try:
                    x = int(linelist[n])
                    isInt = True
                except ValueError:
                    x = float(linelist[n])
                current_list.append(x)
            matrix.append(current_list)

    # We calculate the mean value of the diagonal entries.
    mean = getMean(matrix, isInt)
    output_open = open(output, "w")

    # We write the new matrix to the output.
    for elem in amino_acids:
        output_open.write(elem + " ")

    output_open.write("\n")

    for n in range(len(matrix)):
        output_open.write(amino_acids[n] + " ")
        for m in range(len(matrix[n])):
            if n == m:
                output_open.write(str(mean) + " ")
            else:
                output_open.write(str(matrix[n][m]) + " ")
        output_open.write("\n")


if __name__ == '__main__':
    main()
