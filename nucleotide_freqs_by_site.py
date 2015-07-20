import os
from optparse import OptionParser
import numpy as np

'''
Usage
-----
python nucleotide_freqs_by_site.py -f FASTADIR -o OUTDIR [-c CODE -g GAPCODE
                                    -d -t DIVOUT]
-f --fastadir
    Path to directory containing fasta files

-o --outdir
    Path to directory in which to store output files of base counts
    for each FASTA file.
    Default = ./out

-c --code 
    String containing the base codes to count in the FASTA
    file, excluding the code for gaps.
    Default = "ACGT".

-g --gapcode
    String containing the symbol used to represent gaps
    in the aligned FASTA file.
    Default = "-". 

-d --calculate-diversity
    If -d is included, calculate a nucleotide diversity estimate for each file.
   
-t --divout
    Output file for diversity estimates.  Ignored if -d is not selected.
    Default = "diversity.tsv"


Input
-----
Directory containing aligned FASTA files

Output
-----
Directory of tables providing the number of occurances of each base at each
position in the alignment.

If -d is specified, a single table is also produced summarising the nucleotide
site diversity for each fasta file.
Nucleotide site divesity is calculated as:
    div = X/(((D*(D-1))/2)*L)
    where
    X is the total number of bases over all positions in the aligned
    FASTA which represent a minor allele (bases which are not the most common
    at a certain position).
    D is the number of sequences in the FASTA file
    L is the length of the alignment in base pairs

'''


def make_fastagrid(fasta):
    '''
    Converts the FASTA file to a numpy array where each row corresponds to
    a sequence and each column to a position.
    '''
    seq = []
    x = 0
    with open(fasta) as input:
        for line in input:
            line = line.strip()
            if x == 0:
                x += 1
            elif line[0] == ">":
                seq = "".join(seq).upper()
                npseq = np.array(list(seq))
                if x == 1:
                    L = len(seq)
                    seqgrid = npseq
                    x += 1
                else:
                    assert len(seq) == L, (
                        "sequence lengths in %s  are not equal" % fasta)
                    seqgrid = np.vstack([seqgrid, npseq])
                seq = []
            else:
                seq.append(line)
    seq = "".join(seq).upper()
    seqgrid = np.vstack([seqgrid, npseq])
    return seqgrid


def score_fastagrid(grid, code):
    '''
    Counts each base in "code" at each position in the output
    from make_fastagrid.
    '''
    Tgrid = grid.T
    scoregrid = []
    rowno = 1
    for line in Tgrid:
        tot = 0
        scores = [rowno]
        for opt in code:
            x = len(np.where(line == opt)[0])
            tot += x
            scores.append(x)
        scores.append(len(line) - tot)
        scoregrid.append(scores)
        rowno += 1
    scoregrid = np.vstack(scoregrid)
    return scoregrid


def calculate_X(scoregrid):
    '''
    Calculates X - the total number of bases over all positions in the aligned
    FASTA which represent a minor allele (bases which are not the most common
    at a certain position).
    Also returns L - alignment length and D - alignment depth
    '''
    X = 0
    L = len(scoregrid) - 1
    D = sum(scoregrid[1, 1:])
    subgrid = scoregrid[:, 1:-2]
    for line in subgrid:
        x = sum(line[np.where(line != max(line))[0]])
        X += x
    return (float(X), float(L), float(D))


def run_all(options):
    '''
    Imports options, runs make_fastagrid, score_fastagrid, calculate_X, writes
    the output files.
    '''
    fastadir = options.fastadir
    assert fastadir is not None, (
        "Provide a directory of FASTA files using -f")
    fastas = os.listdir(fastadir)

    outdir = options.outdir
    code = list(options.code.upper()) + [options.gapcode]
    gapcode = options.gapcode
    div = options.div
    divout = options.divout

    print "fastadir = %s" % (fastadir)
    print "outdir = %s" % (outdir)
    print "code = %s" % (code)
    print "gapcode = %s" % (gapcode)
    print "calculate-diversity = %s" % (div)
    if div is True:
        print "divout = %s" % (divout)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if div is True:
        do = open(divout, "w")
        do.write("File\tX\tLength\tDepth\tPi\n")

    for fasta in fastas:
        grid = make_fastagrid("%s/%s" % (fastadir, fasta))
        scoregrid = score_fastagrid(grid, code)

        if div is True:
            X, L, D = calculate_X(scoregrid)
            pi = X/(((D*(D-1))/2)*L)
            do.write("%s\t%i\t%i\t%i\t%f\n" % (fasta, X, L, D, pi))

        outfile = "%s/%s.txt" % (outdir, fasta)
        H = "Position\t%s\tOther" % "\t".join(code)
        np.savetxt(outfile, scoregrid, fmt="%i", delimiter="\t",
                   header=H, comments="")
    if div is True:
        do.close()


def main():
    parser = OptionParser()
    parser.add_option("-f", "--fastadir", action="store",
                      type="string", dest="fastadir")
    parser.add_option("-o", "--outdir", action="store",
                      type="string", dest="outdir", default="out")
    parser.add_option("-c", "--code", action="store",
                      type="string", dest="code", default="ACGT")
    parser.add_option("-g", "--gapcode", action="store",
                      type="string", dest="gapcode", default="-")
    parser.add_option("-d", "--calculate-diversity", action="store_true",
                      dest="div")
    parser.add_option("-t", "--divout", action="store",
                      type="string", dest="divout", default="diversity.tsv")    

    (options, args) = parser.parse_args()
    run_all(options)


if __name__ == "__main__":
    main()
