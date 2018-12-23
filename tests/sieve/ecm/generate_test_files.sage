import os

current_sage_file = os.path.join (os.getcwd(), sys.argv[0])
sage_file_dir = os.path.dirname (current_sage_file)

# load cofac_utils.sage from CADO_NFS_SOURCE/sieve/ecm
load_dir = os.path.join (sage_file_dir, "..", "..", "..", "sieve", "ecm")
load_attach_path (os.path.normpath (load_dir))
load cofac_utils.sage

################################################################################

def factor_test_line (outfile, pdata, qdata, pcomment="", qcomment=""):
  if qdata[0] < pdata[0]: # Rename so that p refers to the smaller prime
    qdata, pdata, qcomment, pcomment = pdata, qdata, pcomment, qcomment
  p = pdata[0]
  q = qdata[0]
  N = p*q
  pcomment += " (order=%d, lpf=%d)" % (pdata[2], pdata[3])
  qcomment += " (order=%d, lpf=%d)" % (qdata[2], qdata[3])
  outfile.write ("%d %d %d # %s, %s\n" % (N, p, q, pcomment, qcomment))


def write_factor_test_file (outfile, B1, B2, method, param, minq=40, minp=10000):
  header = "# Created with: write_factor_test_file (%d, %d, %d, %d, %d, %d)\n"
  outfile.write (header % (B1, B2, method, param, minq, minp))

  # Composite number < 2^32
  # A B1-smooth factor # FIXME does not test correctly B1-powersmoothness
  p1 = prime_with_lpf_in_range (minp, minq, B1, method, param)
  # A non-smooth cofactor
  q = prime_with_lpf_in_range (minp, B2+1, 0, method, param)
  factor_test_line (outfile, p1, q, "one B1-smooth factor", "one non-smooth cofactor")

  # A B1, B2-smooth factor, but not B1-smooth
  p2 = prime_with_B1_B2_smooth_order (minp, B1, B2, method, param)
  factor_test_line (outfile, p2, q, "one B1,B2-smooth factor", "one non-smooth cofactor")

  # A B1 and a B1,B2-smooth factor
  factor_test_line (outfile, p1, p2, "one B1-smooth factor", "one B1,B2-smooth factor")

  # Find two B1-smooth factors with different power of 2 in the order
  # Backtracking does not work reliably for ECM as in the addition chain with a
  # point of small order, an addition may be used incorrectly where a doubling
  # would be required, causing a zero coordinate before a backtracking
  # checkpoint is reached.
  if method < 2: # ie P-1 or P+2
    p3 = p1
    while p1[2].valuation(2) == p3[2].valuation(2):
      p3 = prime_with_lpf_in_range (p3[0]+1, minq, B1, method, param);
    factor_test_line (outfile, p1, p3, " # Two B1-smooth factors with different power of 2 in the order", "")

  ## XXX remove this test
  q1 = prime_with_lpf_in_range (minp, B1, B1+50, method, param)
  q2 = prime_with_lpf_in_range (minp, q1[3] + 500, B2, method, param)
  factor_test_line (outfile, q1, q2, " # Two B1,B2-smooth factors with LPF in differnet giant-steps", "")

  for v in [33, 49, 65, 97, 127, 128, 200]:
    # A non-smooth cofactor such that the product is >2^v[i]
    q = prime_with_lpf_in_range (floor(2^v / p1[0]), B2+1, 0, method, param)
    factor_test_line (outfile, p1, q, "one B1-smooth factor", "one non-smooth cofactor")
    q = prime_with_lpf_in_range (floor(2^v / p2[0]), B2+1, 0, method, param)
    factor_test_line (outfile, p2, q, "one B1,B2-smooth factor", "one non-smooth cofactor")

################################################################################

def write_order_test_file (outfile, pmin, pmax, method, param):
  p = next_prime (pmin-1)
  while p <= pmax:
    try:
      info = get_order_from_method (method, param, p)
      outfile.write ("%d %d %d\n" % (info[0], info[0], info[2]))
    except (ZeroDivisionError, ArithmeticError):
      pass
    p = next_prime (p)

################################################################################

print "### Output directory:", sage_file_dir

##### factor test
B1 = 100
B2 = 1000
methods = {
    (0, "pm1"): [2],
    (1, "pp1"): [ 2/7, 6/5 ],
    (2, "ecm"): [ 10, 11 ],
    (3, "ecmm12"): [ 2, 4 ],
    (4, "ecmm16"): [ 1 ],
    (5, "ecmem12"): [ 2 ]
  }
for m, params_list in methods.items():
  method, m_str = m
  for param in params_list:
    if param in ZZ:
      param_str = str(param)
    elif param in QQ:
      param_str = "%d%d" % (param.numerator(), param.denominator())
    else:
      raise TypeError ("param must be in ZZ or QQ")

    filename = "test_factor_%s_%s_%d_%d.inp" % (m_str, param_str, B1, B2)
    filepath = os.path.join (sage_file_dir, filename)
    print "### Writing %s ... " % filename
    with open (filepath, 'w') as f:
      write_factor_test_file (f, B1, B2, method, param)
    print "### Done"


##### order test
methods = {
    (2, "ecm"): [ 10, 11, 12 ],
    (3, "ecmm12"): [ 4, 5, 6],
    (4, "ecmm16"): [ 1 ],
    (5, "ecmem12"): [ 2, 3, 4 ],
  }
for pmin in [ 1000, 1000000, 1000000000 ]:
  pmax = pmin + 1000
  for m, params_list in methods.items():
    method, m_str = m
    for param in params_list:
      filename = "test_order_%s_%d_%d_%d.inp" % (m_str, param, pmin, pmax)
      filepath = os.path.join (sage_file_dir, filename)
      print "### Writing %s ... " % filename
      with open (filepath, 'w') as f:
        write_order_test_file (f, pmin, pmax, method, param)
      print "### Done"
