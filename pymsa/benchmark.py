import click

from pymsa import PercentageOfNonGaps, PercentageOfTotallyConservedColumns, read_fasta_file_as_list_of_pairs, \
    print_alignment, SumOfPairs, Entropy, Star, PAM250, Blosum62


@click.command()
@click.option(
    '--input_fasta',
    '-i',
    type=click.File('r'),
    help='Input file in FASTA format.',
)
@click.option(
    '--show',
    '-s',
    default=True,
    help='Whether you want to print the input alignment or not.'
)
def benchmark(input_fasta, show: bool):
    msa = read_fasta_file_as_list_of_pairs(input_fasta.name)

    aligned_sequences = list(pair[1] for pair in msa)
    sequences_id = list(pair[0] for pair in msa)

    if show:
        print_alignment(aligned_sequences, sequences_id)

    # Percentage of non-gaps and totally conserved columns
    non_gaps = PercentageOfNonGaps()
    totally_conserved_columns = PercentageOfTotallyConservedColumns()

    percentage = non_gaps.compute(aligned_sequences)
    print("Percentage of non-gaps: {0} %".format(percentage))

    conserved = totally_conserved_columns.compute(aligned_sequences)
    print("Percentage of totally conserved columns: {0}".format(conserved))

    # Entropy
    value = Entropy().compute(aligned_sequences=aligned_sequences)
    print("Entropy score: {0}".format(value))

    # Sum of pairs
    value = SumOfPairs(Blosum62()).compute(aligned_sequences=aligned_sequences)
    print("Sum of Pairs score (Blosum62): {0}".format(value))

    value = SumOfPairs(PAM250()).compute(aligned_sequences=aligned_sequences)
    print("Sum of Pairs score (PAM250): {0}".format(value))

    # Star
    value = Star(Blosum62()).compute(aligned_sequences=aligned_sequences)
    print("Star score (Blosum62): {0}".format(value))

    value = Star(PAM250()).compute(aligned_sequences=aligned_sequences)
    print("Star score (PAM250): {0}".format(value))


if __name__ == '__main__':
    benchmark()
