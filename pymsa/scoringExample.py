import logging

from pymsa.score import Score, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, \
    SumOfPairs
from pymsa.substitutionmatrix import PAM250, Blosum62

"""
This program is intended to show some examples of using the scores in pyMSAScoring.
"""


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def compute_score(score: Score, sequences: list) -> float:
    """ This function applies an score to a multiple sequence alignment
    :param score: Score method
    :param sequences: List of sequences
    :return: The score value
    """
    return score.compute(sequences)


def run_all_scores(msa: list) -> None:
    sequences = list(pair[1] for pair in msa)

    non_gaps = PercentageOfNonGaps()
    totally_conserved_columns = PercentageOfTotallyConservedColumns()

    percentage = non_gaps.compute(sequences)
    conserved = totally_conserved_columns.compute(sequences)
    logger.info("Percentage of non-gaps: {0} %".format(percentage))
    logger.info("Percentage of totally conserved columns: {0}".format(conserved))

    score_method = Entropy()
    value = compute_score(score_method, sequences)
    logger.info("Entropy score: {0}".format(value))

    score_method = SumOfPairs(Blosum62())
    value = compute_score(score_method, sequences)
    logger.info("SumOfPairs score (Blosum62): {0}".format(value))

    score_method = SumOfPairs(PAM250())
    value = compute_score(score_method, sequences)
    logger.info("{0} score (PAM250): {1}".format(score_method.get_name(), value))

    score_method = Star(Blosum62())
    value = compute_score(score_method, sequences)
    logger.info("Star score (Blosum62): {0}".format(value))

    score_method = Star(PAM250())
    value = compute_score(score_method, sequences)
    logger.info("Star score (PAM250): {0}".format(value))

if __name__ == '__main__':
    msa = [("s1", "ACTG"), ("S2", "A-T-")]
    run_all_scores(msa)
