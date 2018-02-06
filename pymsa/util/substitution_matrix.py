import re


def read_matrix_from_file(path_to_file: str) -> dict:
    distance_matrix = {}
    header = ()

    try:
        with open(path_to_file, 'r') as matrix:
            for line in matrix.readlines():
                if not line.startswith('#') and line.strip():
                    # remove leading and trailing spaces and then replace consecutive whitespace characters
                    tmp = re.sub("\s+", " ", line.strip()).split(" ")

                    # if not header was specified, use the first line from the input file
                    if line.startswith(' '):
                        header = tmp
                    else:
                        for i in range(len(header) - 1):
                            if (tmp[0], header[i + 1]) not in matrix and (header[i + 1], tmp[0]) not in matrix:
                                distance_matrix.update({(tmp[0], header[i]): int(tmp[i + 1])})
    except FileNotFoundError:
        raise Exception('File {0} not found!'.format(path_to_file))

    return distance_matrix
