# Usage:  polsel_params("../parameters/factor/")

import sys,os,re,glob
re_param_file = 'params.c*'
default_nq = 1000
default_incr = 60
old_nadall = 0
old_P = 0
old_time = 0
old_n = 0
time_base=2000000
time_p=1.2 # P-collision for each P=2e6 take about 1.2s
time_q=0.02  #  q-collision for each P=2e6 take about 0.02s
"""
On a single core of Intel Q9400 @ 2.66GHz, for RSA155, the
following timings are obtained.

-------------------------------------------
P       |   Each P   |   Each q
-------------------------------------------
5e3     |    nil     |   nil
1e4     |    4ms     |   nil
1e5     |    50ms    |   nil
2e5     |    130ms   |   nil
5e5     |    316ms   |   4ms
1e6     |    628ms   |   8ms
2e6     |    1244ms  |   24ms
4e6     |    2472ms  |   56ms
8e6     |    4900ms  |   120ms
1e7     |    6130ms  |   144ms
--------------------------------------------

For larger N, each P may take more time, while each q could be similar.
We use p2e6 as a benchmark. Note this doesn't consider ropt time.
"""

def polsel_params(fdpath):
    global old_nadall
    old_nadall = 0
    """ feed a folder """

    files = glob.glob(os.path.join(fdpath, re_param_file))
    polsel_params_sort(files)

    for fpath in files:
        # print "# ", fpath
        polsel_param(fpath)
    return

def polsel_param(fpath):
    """ feed a single file
    degree=4
    polsel_P=7400
    polsel_maxnorm=36.0
    polsel_admax=1e6
    polsel_incr=210
    polsel_lq=3
    """
    global default_incr
    flag = 0
    incr = default_incr
    fi = open(fpath,"r")
    s = fi.readlines()
    admin = 0

    for i in range(len(s)):
        Hdegree = re.match("tasks.polyselect.degree[ ]*=[ ]*(\d)", s[i])
        Hnq = re.match("tasks.polyselect.nq[ ]*=[ ]*(\d+)", s[i])
        Hp = re.match("tasks.polyselect.P[ ]*=[ ]*(\d+)", s[i])
        Hincr = re.match("tasks.polyselect.incr[ ]*=[ ]*(\d+)", s[i])
        Hadmax = re.match("tasks.polyselect.admax[ ]*=[ ]*([-+]?[0-9]*\.?[0-9]+[eE]?[0-9]+)", s[i])
        Hadmin = re.match("tasks.polyselect.admin[ ]*=[ ]*([-+]?[0-9]*\.?[0-9]+[eE]?[0-9]+)", s[i])
        Hadrange = re.match("tasks.polyselect.adrange[ ]*=[ ]*([-+]?[0-9]*\.?[0-9]+[eE]?[0-9]+)", s[i])

        if Hdegree != None:
            degree = ZZ(Hdegree.groups()[0])
            flag += 1
        if Hp != None:
            p = ZZ(Hp.groups()[0])
            flag += 1
        if Hnq != None:
            nq = ZZ(Hnq.groups()[0])
            flag += 1
        if Hadmax != None:
            admax = RR(Hadmax.groups()[0])
            flag += 1
        if Hadrange != None:
            adrange = RR(Hadrange.groups()[0])
            flag += 1
        if Hadmin != None:
            print Hadmin.groups()[0]
            admin = RR(Hadmin.groups()[0])
        if Hincr != None:
            incr = ZZ(Hincr.groups()[0])
            if (incr==0):
                sys.exit("Error: -incr can't be 0" + fpath)                

        if (flag == 5):
            polsel_param_check (fpath, degree, p, nq, admin, admax, adrange, incr)
            break

    # mandatory flags are degree, p, nq, admax and adrange
    # admin is optional (default 0)
    if (flag < 5):
        sys.exit("Error: some parameters are missing in " + fpath)

    fi.close()
    return

def polsel_param_check(fpath, degree, p, nq, admin, admax, adrange, incr):
    """ check parameters """
    global old_nadall
    global old_P
    global old_n
    global old_time
    global default_nq
    global time_p
    global time_q
    global time_base

    # print info
    print "#", os.path.basename(fpath),
    print "d="+str(degree),
    print "nq="+str(nq),
    print "admin="+str('%.2e'%admin),
    print "admax="+str('%.2e'%admax),
    print "adrange="+str('%.2e'%adrange),
    print "incr="+str(incr),
    print "P="+str(p)

    admax -= admin # since we consider admax-admin below

    # check whether N is too small
    fname = os.path.basename(fpath)
    ndigits = ZZ(re.match("params.c(\d+)", fname).groups()[0])
    if (ndigits < 30):
        sys.exit("Error: number N too small in " + fpath + "\n")

    # check if nad is increasing monotonically
    nadall = RR(admax/incr)*nq
    exp_collisions = nadall/(2*log(p)^2)
    exp_time = (time_p*admax/incr + time_q*(admax/incr)*nq)*p/time_base
    print "# num_ad="+str('%.2e'%nadall),
    print "exp_coll="+str('%.2e'%(exp_collisions)),
    print "exp_time="+str('%.2e'%(exp_time))+"s"
    
    # check if P is increasing monotonically
    if p < old_P:
        print "--> Error: P for last (",
        print old_n, str(old_P),
        print ") is larger than for (",
        print ndigits, str(p), ")" + "\n"
    old_P = p

    # check admax/adrange is at least 100, so that we can split the work into
    # several threads
    if admax / adrange < 100:
        print "--> Warning: admax / adrange < 100: " + str(admax / adrange)

    # check adrange/incr <= 200, otherwise we risk a timeout
    # several threads
    if adrange / incr > 200:
        print "--> Warning: adrange / incr > 200: " + str(adrange / incr)

    if (nadall < old_nadall):
        print "--> Error: effort for last (",
        print old_n, str('%.2e'%old_nadall),
        print ") is larger than for (",
        print ndigits, str('%.2e'%nadall), ")" + "\n"
    else:
        print ""
    old_nadall = nadall
    old_n = ndigits

def tryint(c):
    try:
        return int(c)
    except:
        return c

def polsel_params_sort_method(s):
    return [tryint(c) for c in re.split('([0-9]+)', s)]

def polsel_params_sort(files):
    """ sort file name by size of number """
    files.sort(key=polsel_params_sort_method)
