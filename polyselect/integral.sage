
attach linear_interpolation.sage

def pdf_normal(x):
    """
    Returns the probability density function (pdf) of the
    standard normal distribution.
    """
    return exp(-x^2/2)/sqrt(2*RR.pi())

def cdf_normal(x):
    """
    Returns the cumulative distribution function (cdf) of the
    standard normal distribution.
    """
    return 1-1/2*error_fcn(1.0*float(x/sqrt(2)))


def pdf_normal_minimum(n, x):
    """
    returns the probability density function (pdf) of the first-order
    statistics of the standard normal distribution.
    """
    return n*pdf_normal(x)*(1-cdf_normal(x))^(n-1)

# note: numerical_integral(lambda x: x*pdf_normal_minimum(2^i,x), -infinity,+infinity)
# can be approximated asymptotically by
# -sqrt(2*ln(2)*i)+(ln(ln(2)*i)+ln(4*pi)-2*0.5772)/(2*sqrt(2*ln(2)*i))
# where a relatively good fit of the gap between this expansion and the
# real value for reasonably small terms yields:
# -sqrt(2*ln(2)*i)-(ln(ln(2)*i)+ln(4*pi)-2*0.5772)/(2*sqrt(2*ln(2)*i))-0.1727*exp(-0.2759*ln(2)*i)
# see http://gupea.ub.gu.se/bitstream/2077/3092/1/correction.pdf
# The cited book by Cramer contains more info (page 374).
def build_minimum_normal_table(M,B):
    """
    precomputes a table with the mean and standard deviation of the
    first-order statistics of the standard normal distribution. This
    function takes some time, so it is probably best to leave its result
    precomputed and set aside.
    B indicates the bound integration interval (20 is a good start)
    10^M is the largest sample size to be considered.
    The table is returned with accuracy 0.1 for log(sample_size)
    """
    table=[]
    for p in range(0,M*10):
        B=100
        n=10^(p/10)
        m=numerical_integral(lambda x: (x*pdf_normal_minimum(n,x)),[-B,B])
        s=numerical_integral(lambda x: (x^2*pdf_normal_minimum(n,x)),[-B,B])
        e=sqrt(s[0]-m[0]^2)
        table.append((m[0], e))
    return table

# This table contains the value of expected_alpha(10^(x/10)) for x below 140
# It corresponds to a sieve space below 10^14
# table=build_minimum_normal_table(14,20)
minimum_normal_table=[                                                  \
( 0.000000,1.000000), (-0.200613,0.935269), (-0.387600,0.877603), 	\
(-0.562447,0.826146), (-0.726462,0.780145), (-0.880795,0.738939), 	\
(-1.026459,0.701946), (-1.164348,0.668654), (-1.295245,0.638619), 	\
(-1.419844,0.611451), (-1.538753,0.586808), (-1.652510,0.564395), 	\
(-1.761591,0.543951), (-1.866416,0.525252), (-1.967357,0.508101), 	\
(-2.064744,0.492324), (-2.158871,0.477773), (-2.249998,0.464315), 	\
(-2.338359,0.451836), (-2.424162,0.440235), (-2.507594,0.429424), 	\
(-2.588822,0.419324), (-2.667997,0.409867), (-2.745256,0.400993), 	\
(-2.820721,0.392648), (-2.894505,0.384784), (-2.966707,0.377360), 	\
(-3.037420,0.370337), (-3.106730,0.363683), (-3.174711,0.357367), 	\
(-3.241436,0.351362), (-3.306968,0.345645), (-3.371368,0.340194), 	\
(-3.434690,0.334989), (-3.496986,0.330013), (-3.558302,0.325250), 	\
(-3.618682,0.320684), (-3.678167,0.316304), (-3.736794,0.312096), 	\
(-3.794600,0.308050), (-3.851616,0.304156), (-3.907874,0.300405), 	\
(-3.963402,0.296788), (-4.018228,0.293297), (-4.072378,0.289926), 	\
(-4.125874,0.286667), (-4.178741,0.283514), (-4.230999,0.280463), 	\
(-4.282669,0.277507), (-4.333770,0.274642), (-4.384319,0.271863), 	\
(-4.434335,0.269165), (-4.483833,0.266546), (-4.532829,0.264001), 	\
(-4.581337,0.261527), (-4.629373,0.259121), (-4.676948,0.256779), 	\
(-4.724076,0.254500), (-4.770769,0.252279), (-4.817040,0.250115), 	\
(-4.862897,0.248005), (-4.908354,0.245947), (-4.953419,0.243939), 	\
(-4.998103,0.241979), (-5.042415,0.240065), (-5.086364,0.238195), 	\
(-5.129958,0.236368), (-5.173206,0.234582), (-5.216117,0.232835), 	\
(-5.258697,0.231127), (-5.300954,0.229455), (-5.342895,0.227819), 	\
(-5.384528,0.226216), (-5.425858,0.224647), (-5.466892,0.223110), 	\
(-5.507637,0.221604), (-5.548098,0.220127), (-5.588282,0.218679), 	\
(-5.628192,0.217260), (-5.667836,0.215867), (-5.707218,0.214500), 	\
(-5.746344,0.213159), (-5.785218,0.211843), (-5.823845,0.210550), 	\
(-5.862229,0.209281), (-5.900375,0.208034), (-5.938288,0.206809), 	\
(-5.975971,0.205605), (-6.013429,0.204422), (-6.050666,0.203258), 	\
(-6.087685,0.202115), (-6.124490,0.200990), (-6.161084,0.199884), 	\
(-6.197472,0.198795), (-6.233657,0.197725), (-6.269642,0.196671), 	\
(-6.305430,0.195633), (-6.341024,0.194612), (-6.376428,0.193606), 	\
(-6.411645,0.192616), (-6.446677,0.191638), (-6.481527,0.190678), 	\
(-6.516198,0.189727), (-6.550694,0.188791), (-6.585015,0.187864), 	\
(-6.619166,0.186953), (-6.653148,0.186047), (-6.686964,0.185168), 	\
(-6.720613,0.184352), (-6.754103,0.183502), (-6.787434,0.182662), 	\
(-6.820609,0.181827), (-6.853628,0.181008), (-6.886504,0.180044), 	\
(-6.919215,0.179334), (-6.951789,0.178404), (-6.984211,0.177681), 	\
(-7.016493,0.176777), (-7.048618,0.176200), (-7.080613,0.175395), 	\
(-7.112464,0.174768), (-7.144195,0.173778), (-7.175781,0.172935), 	\
(-7.207227,0.172275), (-7.238548,0.171421), (-7.269712,0.171191), 	\
(-7.300762,0.170661), (-7.331666,0.170531), (-7.362484,0.169518), 	\
(-7.393083,0.170567), (-7.423669,0.169378), (-7.454087,0.168451), 	\
(-7.484383,0.169029), (-7.514552,0.169017), (-7.544737,0.166001), 	\
(-7.574306,0.174355), (-7.604172,0.173667), (-7.634136,0.169100), 	\
(-7.663866,0.166641), (-7.693547,0.160234), (-7.722988,0.159124) ]

def build_reduced_minimum_normal_table(mean,sdev,epsilon):
    rt0=reduced_interpolation_table(minimum_normal_table,epsilon/sdev)
    rt=[(z[0],z[1]*sdev+vector([mean,0])) for z in rt0]
    return rt

# Since minimum_normal_table has indexes scaled by 10, these tables as
# well.
# reduced_minimum_normal_table=build_reduced_minimum_normal_table(0,1,0.05)
# reduced_alpha_table=build_reduced_minimum_normal_table(0,0.95,0.05)

def print_c(table):
    """
    Prints the table for inclusion in a C program.
    """
    st=["{ %d, %f, %f },\n" % (z[0],z[1][0],z[1][1]) for z in table]
    print (reduce(operator.add, st))


# print create_c_table(table,0.005,0.95)
# { 0, 0.000000, 0.950000 },
# { 1, -0.190582, 0.888506 },
# { 2, -0.368220, 0.833723 },
# { 3, -0.534325, 0.784839 },
# { 5, -0.836755, 0.701992 },
# { 7, -1.106131, 0.635221 },
# { 9, -1.348852, 0.580878 },
# { 11, -1.569884, 0.536175 },
# { 13, -1.773095, 0.498989 },
# { 15, -1.961507, 0.467708 },
# { 18, -2.221441, 0.429244 },
# { 21, -2.459381, 0.398358 },
# { 24, -2.679685, 0.373016 },
# { 27, -2.885549, 0.351820 },
# { 30, -3.079364, 0.333794 },
# { 34, -3.322137, 0.313512 },
# { 38, -3.549954, 0.296491 },
# { 42, -3.765232, 0.281949 },
# { 46, -3.969804, 0.269338 },
# { 51, -4.212618, 0.255707 },
# { 56, -4.443101, 0.243940 },
# { 61, -4.662936, 0.233650 },
# { 72, -5.115302, 0.214905 },
# { 84, -5.569118, 0.198817 },
# { 97, -6.023973, 0.184881 },
# { 110, -6.448062, 0.173529 },
# { 124, -6.876621, 0.162850 },
# { 140, -7.337067, 0.146136 },

# This evaluates a function as given by its piecewise linear
# approximation
def eval_dichotomy(t,x):
    a=0
    b=len(t)-1
    assert(t[a][0] <= x)
    if x >= t[b][0]:
        # Do something ugly.
        b=-1
        a=-2
    else:
        while b-a > 1:
            c=floor((a+b)/2)
            if t[c][0] <= x:
                a=c
            else:
                b=c
    l=(x-t[a][0])/(t[b][0]-t[a][0])
    v=t[a][1]+l*(t[b][1]-t[a][1])
    return v

# This should exhibit no visible difference.
# reduced_minimum_normal_table=build_reduced_minimum_normal_table(0,1,0.05)
# plot([lambda x: quick_eval(minimum_normal_table,10*x)[0], lambda x: eval_dichotomy(reduced_minimum_normal_table,10*x)[0]],0,17)

