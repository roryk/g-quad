import regex
from collections import Counter
from Bio import SeqIO
from argparse import ArgumentParser
import re

HEADER = "\t".join(["id", "pg4s_and_scores", "total_score",
                    "length_normalized_score"])

def get_pg4s(seq):
    pattern = r'(([gG]{3,}[aAcCtTgGuU]{1,7}){3,}[gG]{3,})'
    matches = re.finditer(pattern, seq)
    context = get_context(seq)
    return filter(lambda x: x, map(context, matches))

def get_context(seq):
    # return a function that extracts the genomic context of a match
    def context(match):
        CONTEXT_LENGTH = 50
        context_start = match.span()[0] - CONTEXT_LENGTH
        context_end = match.span()[1] + CONTEXT_LENGTH
        if context_start < 0 or context_end > len(seq):
            return None
        return seq[context_start:context_end]
    return context

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

    with open(args.fasta) as in_handle:
        print HEADER
        for record in SeqIO.parse(in_handle, "fasta"):
            pg4s = get_pg4s(str(record.seq))
            scores = map(score_sequence, pg4s)
            total = sum(scores)
            print "\t".join(map(str, [record.id, zip(pg4s, scores), total,
                                  total / len(record.seq)]))
