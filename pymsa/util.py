def read_fasta_file_as_list_of_pairs(file_name: str, directory: str='./') -> list:
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
