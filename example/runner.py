import logging

from pymsa.core.score import Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, \
    SumOfPairs, Strike
from pymsa.core.substitution_matrix import PAM250, Blosum62, FileMatrix

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
    value = Entropy().compute(align_sequences=align_sequences)
    logger.info("Entropy score: {0}".format(value))

    # Sum of pairs
    value = SumOfPairs(Blosum62()).compute(align_sequences=align_sequences)
    logger.info("SumOfPairs score (Blosum62): {0}".format(value))

    value = SumOfPairs(PAM250()).compute(align_sequences=align_sequences)
    logger.info("SumOfPairs score (PAM250): {0}".format(value))

    value = SumOfPairs(FileMatrix('PAM380.txt')).compute(align_sequences=align_sequences)
    logger.info("SumOfPairs score (PAM380): {0}".format(value))

    # Star
    value = Star(Blosum62()).compute(align_sequences=align_sequences)
    logger.info("Star score (Blosum62): {0}".format(value))

    value = Star(PAM250()).compute(align_sequences=align_sequences)
    logger.info("Star score (PAM250): {0}".format(value))

    # STRIKE
    value = Strike().compute(align_sequences=align_sequences,
                             sequences_id=sequences_id,
                             chains=['A', 'E', 'A', 'A'])
    logger.info("STRIKE score: {0}".format(value))


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
