from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62, FileMatrix
from pymsa.util.fasta import print_alignment


def run_all_scores(sequences: list) -> None:
    aligned_sequences = list(pair[1] for pair in sequences)
    sequences_id = list(pair[0] for pair in sequences)

    msa = MSA(aligned_sequences, sequences_id)
    print_alignment(msa)

    # Percentage of non-gaps and totally conserved columns
    non_gaps = PercentageOfNonGaps(msa)
    totally_conserved_columns = PercentageOfTotallyConservedColumns(msa)

    percentage = non_gaps.compute()
    print("Percentage of non-gaps: {0} %".format(percentage))

    conserved = totally_conserved_columns.compute()
    print("Percentage of totally conserved columns: {0}".format(conserved))

    # Entropy
    value = Entropy(msa).compute()
    print("Entropy score: {0}".format(value))

    # Sum of pairs
    value = SumOfPairs(msa, Blosum62()).compute()
    print("Sum of Pairs score (Blosum62): {0}".format(value))

    value = SumOfPairs(msa, PAM250()).compute()
    print("Sum of Pairs score (PAM250): {0}".format(value))

    value = SumOfPairs(msa, FileMatrix('PAM380.txt')).compute()
    print("Sum of Pairs score (PAM380): {0}".format(value))

    # Star
    value = Star(msa, Blosum62()).compute()
    print("Star score (Blosum62): {0}".format(value))

    value = Star(msa, PAM250()).compute()
    print("Star score (PAM250): {0}".format(value))


if __name__ == '__main__':
    sequences = [("1g41",
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
    run_all_scores(sequences)
