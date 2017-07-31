def read_fasta_file_as_list_of_pairs(directory: str, file_name: str) -> list:
    list_of_pairs = []
    key = ''
    value = ''

    with open(directory+file_name, 'r') as file:
        for line in file:
            if line[0] == '>':
                if key != '':
                    list_of_pairs.append((key, value))
                key = line[1:].rstrip()
                value = ''
            else:
                value += line.rstrip()
    list_of_pairs.append((key, value))
    file.close()

    return list_of_pairs


def save_fasta_serialized(directory: str, file_name: str) -> None:
    with open(directory+file_name+"_serialized.txt", 'w') as output:
        list = read_fasta_file_as_list_of_pairs(directory, file_name)

        for (pair1, pair2) in list:
            output.write(pair1 + ',' + pair2 + '\n')
