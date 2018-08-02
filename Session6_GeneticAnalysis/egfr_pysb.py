from pysb import *
from pysb import pathfinder


# this is your pythonanywhere.com username
user_name = 'rah'

bngl_path = '/home/'+str(user_name)+'/BioNetGen-2.3.1/'
pathfinder.set_path('bng', bngl_path)
pathfinder.get_path('bng')

Model()

# Physical and geometric constants
Parameter('NA', 6.0e23)      # Avogadro's num
Parameter('f', 0.01)         # scaling factor
Expression('Vo', f * 1e-10)  # L
Expression('V', f * 3e-12)   # L

# Initial concentrations
Parameter('EGF_conc', 2e-9)             # nM
Expression('EGF0', EGF_conc * NA * Vo)  # nM
Expression('EGFR0', f * 1.8e5)          # copy per cell

# Rate constants
Expression('kp1', 9.0e7 / (NA * Vo))  # input /M/sec
Parameter('km1', 0.06)                # /sec

Monomer('EGF', ['R'])
Monomer('EGFR', ['L', 'CR1', 'Y1068'], {'Y1068': ['U', 'P']})

Initial(EGF(R=None), EGF0)
Initial(EGFR(L=None, CR1=None, Y1068='U'), EGFR0)

Rule('egf_binds_egfr', EGF(R=None) + EGFR(L=None) | EGF(R=1) % EGFR(L=1), kp1, km1)

# Species LR EGF(R!1).EGFR(L!1)
Observable('Lbound', EGF(R=ANY))  # Molecules

from pysb.integrate import Solver
t = [0, 10, 20, 30, 40, 50, 60]
solver = Solver(model, t)
solver.run()
print(solver.y[:, 1])
import matplotlib.pyplot as plt
plt.plot(solver.y)
# plt.show()
plt.savefig('egfr_observables.png')
plt.close()