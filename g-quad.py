import regex
from collections import Counter
from Bio import SeqIO
from argparse import ArgumentParser
import re


def score_matches(matches):
    score = 0
    for m in matches:
        score += score_sequence(m.group(0))
    return score

def score_sequence(seq):
    """
    implementation of G-quadraplex scoring from this paper
    http://nar.oxfordjournals.org/content/early/2013/10/10/nar.gkt904.full
    """
    G_score = 0
    for i in range(1, len(seq) + 1):
        G_score += len(regex.findall('[gG]{%d}' % i, seq, overlapped=True)) * 10 * i
    C_score = 0
    for i in range(1, len(seq) + 1):
        C_score += len(regex.findall('[cC]{%d}' % i, seq, overlapped=True)) * 10 * i

    if C_score == 0:
        C_score = 1

    return float(G_score) / float(C_score)


if __name__ == "__main__":
    parser = ArgumentParser(description="Calculate the G-quadraplex score for a FASTA file "
            "of sequences. The score is calculated as per this paper: "
            "http://nar.oxfordjournals.org/content/early/2013/10/10/nar.gkt904.full")
    parser.add_argument("fasta", help="A FASTA file of sequences")
    args = parser.parse_args()

    pattern = r'(([gG]{3,}[aAcCtTgGuU]{1,7}){3,}[gG]{3,})'
    with open(args.fasta) as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            matches = re.finditer(pattern, str(record.seq))
            print record.id, score_matches(matches)
