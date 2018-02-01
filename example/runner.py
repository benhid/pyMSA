import logging

from pymsa.core.score import Score, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, \
    SumOfPairs, Strike
from pymsa.core.substitutionmatrix import PAM250, Blosum62

"""
This program is intended to show some examples of using the scores in pyMSA
"""


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def compute_score(score: Score, sequences: list) -> float:
    """ This function applies an core to a multiple sequence alignment
    :param score: Score method
    :param sequences: List of sequences
    :return: The core value
    """
    return score.compute(sequences)


def run_all_scores(msa: list) -> None:
    align_sequences = list(pair[1] for pair in msa)
    sequences_id = list(pair[0] for pair in msa)

    # Percentage of non-gaps and totally conserved columns
    non_gaps = PercentageOfNonGaps()
    totally_conserved_columns = PercentageOfTotallyConservedColumns()

    percentage = non_gaps.compute(align_sequences)
    conserved = totally_conserved_columns.compute(align_sequences)
    logger.info("Percentage of non-gaps: {0} %".format(percentage))
    logger.info("Percentage of totally conserved columns: {0}".format(conserved))

    # Entropy
    score_method = Entropy()
    value = compute_score(score_method, align_sequences)
    logger.info("Entropy core: {0}".format(value))

    # Sum of pairs
    score_method = SumOfPairs(Blosum62())
    value = compute_score(score_method, align_sequences)
    logger.info("SumOfPairs core (Blosum62): {0}".format(value))

    score_method = SumOfPairs(PAM250())
    value = compute_score(score_method, align_sequences)
    logger.info("{0} core (PAM250): {1}".format(score_method.get_name(), value))

    # Star
    score_method = Star(Blosum62())
    value = compute_score(score_method, align_sequences)
    logger.info("Star core (Blosum62): {0}".format(value))

    score_method = Star(PAM250())
    value = compute_score(score_method, align_sequences)
    logger.info("Star core (PAM250): {0}".format(value))

    # STRIKE
    value = Strike().compute(align_sequences, sequences_id, ['A', 'E', 'A', 'A'])
    logger.info("STRIKE core: {0}".format(value))


if __name__ == '__main__':
    msa = [("1g41",
            "S-EMTPREIVSELDQHIIGQADAKRAVAIALRNRWRRMQLQEPLRHE--------VTP-KNILMIGPTGVGKTEIARRLAKLANAPFIKVEATKFT----"
            "VGKEVDSIIRDLTDSAMKLVRQQEIAKNR---------------------------------------------------------------------LI"
            "DDEAAKLINPEELKQKAIDAVE--QNGIVFIDEIDKICKKGEYSGADVSREGVQRDLLPLVEGSTVSTKHGMVKTDHILFIASGAFQVARPSDL------"
            "-----------IPELQGRLPIR-VEL---TALSAADFERILTEPHASLTEQYKALMATEGVNIAFTTDAVKKIAEAAFRVNEKTENIGARRLHTVMERLM"
            "DKISFSASDMNGQTVNIDAAYVADALGEVVENEDLSRFIL"),
           ("1e94",
            "HSEMTPREIVSELDKHIIGQDNAKRSVAIALRNRWRRMQLNEELRHE--------VTP-KNILMIGPTGVGKTEIARRLAKLANAPFIKVEATKFTEVGY"
            "VGKEVDSIIRDLTDAAVKMVRVQAIEKNRYRAEELAEERILDVLIPPAKNNWGQTEQQQEPSAARQAFRKKLREGQLDDKEIEKQKARKLKIKDAMKLLI"
            "EEEAAKLVNPEELKQDAIDAVE--QHGIVFIDEIDKICKRGESSGPDVSREGVQRDLLPLVEGCTVSTKHGMVKTDHILFIASGAFQIAKPSDL------"
            "-----------IPELQGRLPIR-VEL---QALTTSDFERILTEPNASITVQYKALMATEGVNIEFTDSGIKRIAEAAWQVNESTENIGARRLHTVLERLM"
            "EEISYDASDLSGQNITIDADYVSKHLDALVADEDLSRFIL"),
           ("1e32",
            "R-ED-EEESLNEVGYDDVGG--CRKQLAQ-----I-KEMVELPLRHPALFKAIGVKPP-RGILLYGPPGTGKTLIARAVANETGAFFFLINGPEIM-SKL"
            "A-GESESN--------------------------------------------------------------------------------------------"
            "-------------LRKAFEEAEKNAPAIIFIDELDAIAPKREKTHGEVERRIVSQ-LLTLMDGL--------KQRAHVIVMAATN----RPNSIDPALRR"
            "FGRFDREVDIGIPDATGRLEILQIHTKNMKLADDVDLEQVANETHGH---------------------------------------VGADLAALCSEAAL"
            "QAIRKKMDLIDLEDETIDAEVM-NSL-AVTMDDFRWALSQ"),
           ("1d2n",
            "------EDYASYIMNGIIKWGDP---VTRVLD--DGELLVQQTKNSD--------RTPLVSVLLEGPPHSGKTALAAKIAEESNFPFIKICSPDKM-IGF"
            "SETAKCQA--------------------------------------------------------------------------------------------"
            "-------------MKKIFDDAYKSQLSCVVVDDIERLLDYV-PIGPRFSNLVLQA-LLVLLKKA-------PPQGRKLLIIGTTS----R-KDVLQEMEM"
            "LNA---------------------------------FSTTIHVPNIATGEQL--LEALEL-LGNFKDKE---RTTIAQQVKGKKVWIGIKKLLMLIEM--"
            "-------------SLQMDPEYRVRKFLALLREEGAS-PLD")]

    run_all_scores(msa)
