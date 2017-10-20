# Description of the foward finite-difference "algorithm".
from opal.core.algorithm import Algorithm
from opal.core.parameter import Parameter
from opal.core.measure   import Measure

# Define Algorithm object.
LAS = Algorithm(name='LAS', description='Lattice Siever')

# Register executable for LAS.
LAS.set_executable_command('python las_run.py')

# Register parameter file used by black-box solver to communicate with LAS.
#LAS.set_parameter_file('las.param')  # Should be chosen automatically and hidden.

# Define parameter and register it with algorithm.
lim0 = Parameter(kind='integer', default=lim0_def, bound=(lim0_min, lim0_max), name='lim0', description='Factor base bound, side 0')
lim1 = Parameter(kind='integer', default=lim1_def, bound=(lim1_min, lim1_max), name='lim1', description='Factor base bound, side 1')
lpb0 = Parameter(kind='integer', default=lpb0_def, bound=(lpb0_min, lpb0_max), name='lpb0', description='Large prime bound, side 0')
lpb1 = Parameter(kind='integer', default=lpb1_def, bound=(lpb1_min, lpb1_max), name='lpb1', description='Large prime bound, side 1')
mfb0 = Parameter(kind='integer', default=mfb0_def, bound=(mfb0_min, mfb0_max), name='mfb0', description='Cofactorization bound, side 0')
mfb1 = Parameter(kind='integer', default=mfb1_def, bound=(mfb1_min, mfb1_max), name='mfb1', description='Cofactorization bound, side 1')
ncurves0 = Parameter(kind='integer', default=ncurves0_def, bound=(ncurves0_min,ncurves0_max), name='ncurves0', description='Cofactorization curves, side 0')
ncurves1 = Parameter(kind='integer', default=ncurves1_def, bound=(ncurves1_min,ncurves1_max), name='ncurves1', description='Cofactorization curves, side 1')
I = Parameter(kind='integer', default=I_def, bound=(I_min, I_max), name='I', description='Sieve region size')
qmin = Parameter(kind='integer', default=qmin_def, bound=(qmin_min, qmin_max), name='qmin', description='Special-q lower bound')
bkthresh1 = Parameter(kind='integer', default=bkthresh1_def, bound=(bkthresh1_min, bkthresh1_max), name='bkthresh1', description='Level-2 bucket sieve bound')
LAS.add_param(lim0)
LAS.add_param(lim1)
LAS.add_param(lpb0)
LAS.add_param(lpb1)
LAS.add_param(mfb0)
LAS.add_param(mfb1)
LAS.add_param(ncurves0)
LAS.add_param(ncurves1)
LAS.add_param(I)
LAS.add_param(qmin)
LAS.add_param(bkthresh1)

# Define relevant measure and register with algorithm.
sievetime = Measure(kind='real', name='SIEVETIME', description='Time in the sieving')
rels = Measure(kind='integer', name='RELATIONS', description='Relations found in the sieving')
LAS.add_measure(sievetime)
LAS.add_measure(rels)
