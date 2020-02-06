from pymsa.core.msa import MSA


def read_fasta_file_as_list_of_pairs(file_name: str) -> list:
    """
    Read a file in FASTA format as list of pairs (sequence id, sequence).

    :param file_name: FASTA file.
    :return: List of pairs.
    """
    list_of_pairs = []
    key = ''
    value = ''

    with open(file_name, 'r') as file:
        for line in file:
            if line[0] == '>':
                if key != '':
                    list_of_pairs.append((key, value))
                key = line[1:].rstrip()
                value = ''
            else:
                value += line.rstrip()

    list_of_pairs.append((key, value))
    return list_of_pairs


def print_alignment(msa: MSA, cx_point: int = 100):
    sub_sequences = [[]] * msa.number_of_sequences
    for i, sequence in enumerate(msa.sequences):
        sub_sequences[i] = [sequence[i: i + cx_point] for i in range(0, len(sequence), cx_point)]

    for k in range(len(sub_sequences[0])):
        sequences = [item[k] for item in sub_sequences]
        colour_scheme = [0] * len(sequences[0])

        for i, column in enumerate(zip(*sequences)):
            if len(set(column)) <= 1:
                colour_scheme[i] = 1 if len(set(column)) <= 1 else 0
            else:
                if set(column).issubset(['I', 'L', 'V']):
                    colour_scheme[i] = 2
                elif set(column).issubset(['F', 'W', 'Y']):
                    colour_scheme[i] = 2
                elif set(column).issubset(['K', 'R', 'H']):
                    colour_scheme[i] = 2
                elif set(column).issubset(['D', 'E']):
                    colour_scheme[i] = 2
                elif set(column).issubset(['G', 'A', 'S']):
                    colour_scheme[i] = 2
                elif set(column).issubset(['T', 'N', 'Q', 'M']):
                    colour_scheme[i] = 2

        longest_id = len(max(msa.ids, key=len))

        for sequence, id in zip(sequences, msa.ids):
            print(id + ' ' * (longest_id - len(id)), end='\t', flush=True)
            for i in range(len(colour_scheme)):
                if colour_scheme[i] == 1:
                    print('\x1b[44m\x1b[97m' + sequence[i] + '\033[0m', end='', flush=True)
                elif colour_scheme[i] == 2:
                    print('\x1b[46m\x1b[97m' + sequence[i] + '\033[0m', end='', flush=True)
                else:
                    print(sequence[i], end="", flush=True)
            print()
        print()
