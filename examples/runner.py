import sys
from pymsa import Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, \
    SumOfPairs
from pymsa import PAM250, Blosum62, FileMatrix
from pymsa.util.fasta import read_fasta_file_as_list_of_pairs, print_alignment


def run_all_scores(msa: list) -> None:
    aligned_sequences = list(pair[1] for pair in msa)
    sequences_id = list(pair[0] for pair in msa)

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

    value = SumOfPairs(FileMatrix('PAM380.txt')).compute(aligned_sequences=aligned_sequences)
    print("Sum of Pairs score (PAM380): {0}".format(value))

    # Star
    value = Star(Blosum62()).compute(aligned_sequences=aligned_sequences)
    print("Star score (Blosum62): {0}".format(value))

    value = Star(PAM250()).compute(aligned_sequences=aligned_sequences)
    print("Star score (PAM250): {0}".format(value))

    # STRIKE
    #value = Strike().compute(aligned_sequences=aligned_sequences,
    #                         sequences_id=sequences_id,
    #                         chains=['A', 'E', 'A', 'A'])
    # print("STRIKE score: {0}".format(value))


if __name__ == '__main__':
    if len(sys.argv) != 2:
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
    else:
        msa = read_fasta_file_as_list_of_pairs(sys.argv[1])

    run_all_scores(msa)
