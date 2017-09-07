import os
import sys
import re
import time
import itertools

# This code is here for information purposes. It can be used to obtain a
# graphical representation of the area hit by the lattice sieving
# program, as the special-q changes. The main purpose is to compare
# different adjustment strategies for the sieve area (--adjust-strategy
# parameter to las). The sage computer algebra system is used (test done
# with sage 7.5.1). This requires commit 1fa54ed81 of cado-nfs (or newer,
# at least until it gets broken).
#
# Here is an example run for a 155-digit number.
#
# $ ./build/localhost/sieve/las -I 14 -poly /tmp/c155.poly -q0 30950000 -q1 530000000 -lim0 17246818 -lim1 30940618 -lpb0 29 -lpb1 29 -mfb0 62 -mfb1 61 -lambda0 2.24 -lambda1 2.2 -ncurves0 13 -ncurves1 23 -fb /tmp/c155.roots.gz  -t 2  -v --adjust-strategy 0  --random-sample 200 -seed 1 > /tmp/strat0.txt
# $ ./build/localhost/sieve/las -I 14 -poly /tmp/c155.poly -q0 30950000 -q1 530000000 -lim0 17246818 -lim1 30940618 -lpb0 29 -lpb1 29 -mfb0 62 -mfb1 61 -lambda0 2.24 -lambda1 2.2 -ncurves0 13 -ncurves1 23 -fb /tmp/c155.roots.gz  -t 2  -v --adjust-strategy 2  --random-sample 200 -seed 1 > /tmp/strat2.txt
# sage display-lattice-adjust.sage /tmp/c155.poly /tmp/strat0.txt /tmp/strat2.txt

# An arbitrary number of relation files may be put on the command line.
# The code will print comparative results for the different adjustment
# strategies, for the special-q's found in the relation files which were
# sieved with both strategies (other special-q's will produce no output
# at all).

# The pdf files produced can conveniently be merged into unique pdf files
# with bash code such as:
#
# for b in {26..29}  ; do pdfjoin `ls /tmp/q$b-*-*.pdf | sort -t- -k+2n` -o /tmp/c155-q$b.pdf ; done
# pdfjoin  /tmp/c155-q??.pdf -o /tmp/c155-q.pdf

sqside = None
skewness=0
f=[0,0]
lpb=[0,0]
logA=0

# By default we're limiting to 25 q's per bitsize. Put something else, or
# "None" if desired.
maxq_per_bitsize = 25

ZP.<x>=PolynomialRing(Integers())

def parse_poly_file(filename):
    global skewness
    with open(filename, "r") as file:
        for line in file:
            for foo in re.finditer("^skew: ([\d\.]+)", line):
                skewness=float(foo.groups()[0])
            for foo in re.finditer("^# f\(x\) = (.*)$", line):
                f[1] = ZP(foo.groups()[0])
            for foo in re.finditer("^# g\(x\) = (.*)$", line):
                f[0] = ZP(foo.groups()[0])

def homo_eval(f,x,y):
    return sum([f[i]*x^i*y^(f.degree()-i) for i in range(0,f.degree()+1)])

def D(f,qb,lp,x,y):
    v=abs(homo_eval(f,x,y))
    vbits=log(v,2)-qb
    return dickman_rho(vbits/lp)

def graphics_common(f, qbits, lpb, logA):
    global sqside
    print "f=",f
    print "qbits=%d on side %d" % (qbits,sqside)
    print "lpb=",lpb
    print "logA=",logA

    _=var('x y');

    qb=[0,0]
    qb[sqside]=qbits
    def p(x,y):
        return D(f[0],qb[0],lpb[0],x,y)*D(f[1],qb[1],lpb[1],x,y)

    K=logA//2 + ceil(qbits/2) + 1

    C=contour_plot(log(p(x*sqrt(skewness),y/sqrt(skewness)),2),
            (x,-2^K,2^K),
            (y,-2^K,2^K),
            plot_points=250,
            fill=False,
            cmap='hsv',
            contours=25)
    plot_roots=[sum([plot(skewness*x/r,(x,-2^K,2^K),color='black') for r in p.real_roots()]) for p in f]
    c=sum(plot_roots) + C
    c[0].set_zorder(1)
    c.set_axes_range(xmin=-2^K,xmax=2^K,ymin=-2^K,ymax=2^K)
    return (K, c)

def graphics_save(K,G,q):
    imagename="/tmp/q%d.pdf"%prod(q)
    G.save(imagename,xmin=-2^K,xmax=2^K,ymin=-2^K,ymax=2^K)
    print "# saved to %s" % imagename

def process_strategy0(file):
    global sqside
    Q0=matrix(2,2)
    Jx=0
    plist0=[]
    for line in file:
        for foo in re.finditer("^# Sieving side-(\d) q=\d+; rho=\d+; a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); b1=(-?\d+); J=(\d+);$", line):
            side=int(foo.groups()[0])
            assert sqside is None or sqside == side
            sqside = side
            Q0=matrix(2,2,[int(x) for x in foo.groups()[1:5]])
            Jx=int(foo.groups()[5])
        for foo in re.finditer("^(-?\d+),(-?\d+)",line):
            a,b=foo.groups()
            plist0.append((int(a),int(b)))
        for foo in re.finditer("^# (\d+) relation", line):
            nrels = int(foo.groups()[0])
            assert nrels == len(plist0)
            return Q0, Jx, plist0

def process_strategy2(file):
    global sqside
    Q1=matrix(2,2)
    Q0=matrix(2,2)
    logI=0
    plist2=[]
    S=matrix(2,2)
    for line in file:
        for foo in re.finditer("# Initial q-lattice: a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); b1=(-?\d+);$", line):
            Q0 = matrix(2,2,[int(x) for x in foo.groups()])
        for foo in re.finditer("# Adjusting by \[(-?\d+),(-?\d+),(-?\d+),(-?\d+)\], logI=(\d+)", line):
            S=matrix(2,2,[int(x) for x in foo.groups()[0:4]])
            logI=int(foo.groups()[4])
        for foo in re.finditer("# Sieving side-(\d) q=\d+; rho=\d+; a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); b1=(-?\d+);", line):
            side=int(foo.groups()[0])
            assert sqside is None or sqside == side
            sqside = side
            Q1=matrix(2,2,[int(x) for x in foo.groups()[1:5]])
            assert Q1 == S * Q0
        for foo in re.finditer("^(-?\d+),(-?\d+)",line):
            a,b=foo.groups()
            plist2.append((int(a),int(b)))
        for foo in re.finditer("^# (\d+) relation", line):
            nrels = int(foo.groups()[0])
            assert nrels == len(plist2)
            return Q1, logI, plist2

def make_rectangles(qr, tup0, tup2):
    q = int(qr)
    r = int((qr-q)*q)
    D = diagonal_matrix([1/sqrt(skewness),sqrt(skewness)])
    def latcorners(Q, I, J):
        return [vector(v)*Q*D for v in [(-I//2,-J),(-I//2,J),(I//2,J),(I//2,-J)]]

    # We rely on a modified binary which processes adjusting
    # strategies 0 and 2 in turn for each special-q. Since the print flow is not
    # exactly the same for the two strategies, we'll process them one
    # after another.

    # print "Now processing %s" % q

    # First strategy 0.
    logI = (logA+1)//2
    Q0, Jx, plist0 = tup0

    # Find I,J for base rectangle
    I=2^(logI)
    J=2^(logA-logI)
    R0=polygon(latcorners(Q0,I,J),fill=False,color='pink',zorder=2)

    print Jx

    R0 += polygon(latcorners(Q0,I,Jx),fill=False,color='red',zorder=3)

    # We don't add the points yet, because we'll have different colors
    # for the intersection.

    # Now strategy 2
    Q1, logI, plist2 = tup2
    # And now for the adjusted one.
    I=2^logI
    J=2^(logA-logI)
    R2 = polygon(latcorners(Q1,I,J),fill=False,color='blue',legend_label=str(q),zorder=3)

    Pc = set(plist0).intersection(set(plist2))
    # P0 = set(plist0) - set(plist2)
    # P2 = set(plist2) - set(plist0)

    legend="+%d-%d=%d" % (
                len(set(plist2) - set(plist0)),
                len(set(plist0) - set(plist2)),
                len(plist2)-len(plist0))

    print "%d %d %s" % (q, r, legend)

    Pc=[(vector(u)*D) for u in Pc]
    P0=[(vector(u)*D) for u in plist0]
    P2=[(vector(u)*D) for u in plist2]
    Pc+=[-u for u in Pc]
    P0+=[-u for u in P0]
    P2+=[-u for u in P2]
    Pc=[tuple(u) for u in Pc]
    P0=[tuple(u) for u in P0]
    P2=[tuple(u) for u in P2]
    R = R0 + R2

    xP0 = points(P0, color='pink', legend_label=str(len(P0)), zorder=4)
    xP2 = points(P2, color='blue', legend_label="%d [%s]" % (len(P2), legend), zorder=4)
    xPc = points(Pc, color='red', legend_label=str(len(Pc)), zorder=5)
    R += xPc + xP0 + xP2

    return R, legend

def process_one_file(file, dict0, dict2):
    global f
    global lpb
    global logA
    qr=None
    for line in file:
        for foo in re.finditer("^.* -A (\d+)", line):
            if logA == 0:
                logA = int(foo.groups()[0])
                print("# Set logA")
            else:
                assert logA == int(foo.groups()[0])
        for foo in re.finditer("^.* -I (\d+)", line):
            if logA == 0:
                logA = 2*int(foo.groups()[0])-1
                print("# Set logA")
            else:
                assert logA == 2*int(foo.groups()[0])-1
        for foo in re.finditer("-lpb(\d) (\d+)", line):
            s,v = foo.groups()
            if lpb[int(s)] == 0:
                lpb[int(s)]=int(v)
                print("# Set lpb%s"%s)
            else:
                assert lpb[int(s)] == int(v)
        for foo in re.finditer("# Making factor base for polynomial g\(x\) = (.*),$", line):
            if f[0] == 0:
                f[0] = ZP(foo.groups()[0])
                print("# Set f0")
            else:
                assert f[0] == ZP(foo.groups()[0])
        for foo in re.finditer(".*f(\d)\(x\) = (.*)$", line):
            s,v = foo.groups()
            if f[int(s)] == 0:
                f[int(s)] = ZP(v)
                print("# Set f%s"%s)
            else:
                assert f[int(s)] == ZP(v)
        # This is for the new binary.
        for foo in re.finditer("# Now sieving side-(\d+) q=(\d+); rho=(\d+)$",line):
            s,v,r=foo.groups()
            q=[1,1]
            q[int(s)]=int(v)
            q[1-int(s)]=1
            # When the next line says "Initial q-lattice" we're doing
            # strategy 2.
            line=file.next()
            qr=prod(q)+int(r)/prod(q)
        for foo in re.finditer("# Initial q-lattice", line):
            assert qr is not None
            file2=itertools.chain([line], file)
            Q1, logI, plist2 = process_strategy2(file2)
            dict2[qr]=(Q1, logI, plist2)
            qr=None
            ## all_rect.append(process_one_specialq(file2, f, lpb, logA, q)
        # This is for the old binary.
        for foo in re.finditer("# Sieving side-(\d+) q=(\d+); rho=(\d+)",line):
            s,v,r=foo.groups()
            q=[1,1]
            q[int(s)]=int(v)
            q[1-int(s)]=1
            qr=prod(q)+int(r)/prod(q)
            # prepend line in the iterable...
            file2=itertools.chain([line], file)
            adjust_strat=0
            Q0, Jx, plist0 = process_strategy0(file2)
            dict0[qr]=(Q0, Jx, plist0)
            ## all_rect.append(process_one_specialq(file, f, lpb, logA, q))
            qr=None


if __name__ == "__main__":
    import sys

    # First try to see whether we have a polynomial file.
    for filename in sys.argv[1:]:
        if re.match(".*\.poly$", filename):
            if skewness != 0:
                raise RuntimeError("several poly filenames given")
            parse_poly_file(filename)
    if skewness == 0:
        raise RuntimeError("no poly filenames given")

    dict0={}
    dict2={}
    for filename in sys.argv[1:]:
        if re.match(".*\.poly$", filename):
            continue
        print("# Analyzing %s" % filename)
        if re.match(".*\.gz$", filename):
            process_one_file(os.popen("zcat " + filename, "r"), dict0, dict2)
        else:
            process_one_file(open(filename, "r"), dict0, dict2)

    common=set(dict0.keys()).intersection(set(dict2.keys()))
    print "%d special-q's in common (among %d for strategy 0 and %d for strategy 2)" % (len(common), len(dict0), len(dict2))

    all_rect=[]
    # now create the rectangles
    for q in sorted(common):
        R,legend = make_rectangles(q, dict0[q], dict2[q])
        all_rect.append((q,R,legend))

    per_bitsize={}
    for qrRl in all_rect:
        qr,R,l = qrRl
        bitsize=ceil(log(qr,2))
        if bitsize not in per_bitsize:
            per_bitsize[bitsize]=[]
        per_bitsize[bitsize].append(qrRl)

    for b in sorted(per_bitsize.keys()):
        n=per_bitsize[b]
        print "%d q's with bit size %d" % (len(n),b)
        K,C=graphics_common(f,b,lpb,logA)
        for qrRl in n[0:maxq_per_bitsize]:
            qr,R,l = qrRl
            G=R+C
            q=int(qr)
            r=int((qr-q)*q)
            imagename="/tmp/q%d-%d-%d.pdf"%(b,q,r)
            G.save(imagename,xmin=-2^K,xmax=2^K,ymin=-2^K,ymax=2^K)
            print "# saved to %s [%s]" % (imagename,l)
