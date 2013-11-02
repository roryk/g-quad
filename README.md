g-quad
======

Calculates the G-quadraplex score for a FASTA file of sequences. This
works by identifying all non-overlapping G-quadraplexes in a sequence
and scoring them as described in:
http://nar.oxfordjournals.org/content/early/2013/10/10/nar.gkt904.full.

A note: 50 nt up and downstream of putative G-quadraplexes are added
to the G-quadraplex before calculating to G-score to place the
G-quadraplex into context. If that runs off the end of the sequence
then the G-quadraplex is discarded.

Usage
=====
python g-quad.py your_sequences.fa
