#!/usr/bin/python
#
# author: Francois Morain (morain@lix.polytechnique.fr)
#
# if index file contains: aa[n] n, then all pairs (a, b) with
# aa[n-1] < a <= aa[n] are in file n.

def usage(s):
    print("./split.py sorted_ab_file dest nfile")
    print("Split the file into subfiles dest0, dest1, ..."
          " and summary in dest/index")


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1 or sys.argv[1] == "help":
        usage(sys.argv[0])
        exit
    inf = open(sys.argv[1], 'r')
    dest = sys.argv[2]
    nf = int(sys.argv[3])
    nl = 0
    while True:
        line = inf.readline()
        if line == "":
            break
        nl += 1
    inf.close()
    inf = open(sys.argv[1], 'r')
    print("There are", nl, "lines that I will split into", nf, "files")
    # mmax = ceil(nl/nf)
    mmax = nl / nf
    if nl % nf != 0:
        mmax += 1
    print("mmax=", mmax)
    index = open(dest + "index", 'w')
    first = True
    n = 0
    in2 = open(dest + str(n), 'w')
    m = 0
    # equal values of a should be in the same file!
    while True:
        line = inf.readline()
        if line == "":
            break
        tmp = line.split()
        if first:
            first = False
            olda = tmp[0]
            in2.write(line)  # at least one item in the first file
            m = m+1
        else:
            if tmp[0] != olda:
                if m >= mmax:
                    in2.close()
                    index.write(olda + " " + str(n) + "\n")
                    print("File", n, "contains", m, "lines")
                    n = n+1
                    in2 = open(dest + str(n), 'w')
                    m = 0
            in2.write(line)
            olda = tmp[0]
            m = m+1
    in2.close()
    index.write(olda + " " + str(n) + "\n")
    print("File", n, "contains", m, "lines")
    index.close()
