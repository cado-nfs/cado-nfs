import sys
import os
# Define a parameter optimization problem in relation to the 
# LAS algorithm.
from las_declaration import LAS

from opal import ModelStructure, ModelData, Model

# set to >1 for no mpi
MPI = int(os.environ.get('NOMAD_MPI', 1))

if MPI > 1:
    from opal.Solvers import NOMADMPI as NOMAD
else:
    from opal.Solvers import NOMAD

# Return the error measure.
def total_time(parameters, measures):
    with open("measures", "a") as f:
        f.write("parameters = %r, measures = %r\n" % (parameters, measures))
    # Now with numerical integration, we want to minimize total time, not the time/rel
    # return sum(map(float, measures["TIME"])) / sum(map(float, measures["RELATIONS"]))
    return sum(map(float, measures["SIEVETIME"]))

# Define parameter optimization problem.
data = ModelData(LAS)
struct = ModelStructure(objective=total_time)  # Unconstrained
model = Model(modelData=data, modelStructure=struct)

# Solve parameter optimization problem.
NOMAD.set_parameter(name='DISPLAY_STATS',
                    value='%3dBBE  %7.1eSOL  %8.3eOBJ  %5.2fTIME %5.2SIEVETIME')
# to limit the number of evaluations, define the environment variable
# NOMAD_MAX_BB_EVAL, for example:
# NOMAD_MAX_BB_EVAL=100 ./optimize.sh ...
import os
max_bb_eval = os.getenv("NOMAD_MAX_BB_EVAL")
if max_bb_eval != None:
   NOMAD.set_parameter(name='MAX_BB_EVAL', value=int(max_bb_eval))
if MPI > 1:
    NOMAD.set_mpi_config("np", MPI)
NOMAD.solve(blackbox=model)
