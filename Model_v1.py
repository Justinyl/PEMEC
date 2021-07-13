import pandas as pd
from pyomo.environ import *

# ================= Initialization =========================

def build_memcap(sl = 20):
    """
    sl = Number of control volumes in the system

    """
    m = ConcreteModel()

    # ========= sets =================
    m.i = Set(initialize=['H2O', 'O2' ,'H2'], doc='set of species')
    m.ci = Set(within = m.i, initialize = ['H2O', 'H2'], doc='subset of species in cathode')
    # creating sl control volumes
    
    m.n = Set(initialize=[i for i in range(0, sl+1)], ordered=True)
  
    # ========= Variables =================
    m.Ca = Var(m.i, m.n, domain = NonNegativeReals, doc='Concentration on the anode side (mol/m^3)')
    m.Cc = Var(m.ci, m.n, domain = NonNegativeReals, doc='Concentration on the cathode side(mol/m^3)')

    m.Nd = Var(m.n, doc = 'diffusion flux of water(mol/m^2*s)')
    m.Neo = Var(m.n, doc = 'electro-osmosis flux of water(mol/m^2*s)')
    m.Nper = Var(m.n, doc = 'Permeation flux of water(mol/m^2*s)') 
    m.Nrxn = Var(m.i, m.n, doc = 'reaction flux(mol/m^2*s)') 

    m.E = Var(m.n, doc = 'open circuit voltage(V)')
    m.dG = Var(m.n, doc = 'Gibbs free energy(J/mol)')
    m.eta = Var(m.n, doc = 'ohmic resistance(ohm)')
    m.j = Var(m.n, domain = NonNegativeReals, doc = 'current density(A/cm^2)')
    # m.V = Var(m.n, domain = NonNegativeReals, doc = 'cell voltage (V)')

    # ========= Auxilary variables =================

    m.ne = Var(m.n, doc = 'Drag coefficient of electro-osmosis(unitless)')
    m.dCadz = Var(m.i, m.n, domain = NonNegativeReals, doc='Partial of concentration wrt z on anode side')
    m.dCcdz = Var(m.ci, m.n, domain = NonNegativeReals, doc='Partial of concentration wrt z on cathode side')

    # ========= Parameters ================
    # Constant 
    m.Tm = Param(initialize = 328, doc = 'Temperature of the membrane') 
    m.L = Param(initialize=0.3, doc='Lenght of the channel (m)') 
    m.Pa0 = Param(initialize=2, doc='inlet 2 (bar)') # Which pressure to use
    m.Pc0 = Param(initialize=100, doc='inlet 2 (bar)') # Which pressure to use
    m.b = Param(initialize=0.105, doc= '0.105 (m)')
    m.u_in = Param(initialize = 3/60, doc = 'velocity(m/s)')     
    m.h_ch = Param(initialize = 0.003, doc = 'height of the channel (m)') 
    
    m.h = Param(initialize = 0.0002, doc = 'membrane thickness (m)')  
    m.F = Param(initialize = 96485, doc = 'Faradays constant (C/mol)')    
    m.gamma = Param(initialize = 0.2, doc = 'reaction fracction of decomposition(unitless)') 
    m.D_H2O = Param(initialize = 1.28e-10, doc = 'diffusivity of H2O(m^2/s)')

    def V_init(m, n):
        return 2
    def j_init(m, n):
        return 1.2
    def ua_init(m, n):
        return value(m.u_in)
    def uc_init(m, n):
        return value(m.u_in)
    def Pa_init(m, n):
        return value(m.Pa0)
    def Pc_init(m, n):
        return value(m.Pc0)
    
    # m.j = Param(m.n, initialize = j_init, doc = 'j_init')
    m.V = Param(m.n, initialize = V_init, doc = 'v_init')
    m.j_avg = Param(initialize = 1, doc = 'Average current density')
    m.ua = Param(m.n, initialize = ua_init, doc = 'Anode velocity initialization' )
    m.uc = Param(m.n, initialize = uc_init, doc = 'Cathode velocity initialization')
    m.Pa = Param(m.n, initialize = Pa_init, doc = 'Anode Pressure (bar) initialization')
    m.Pc = Param(m.n, initialize = Pc_init, doc = 'Cathode Pressure (bar) initialization')

    # Coefficient 
    m.D_H2 = Param(initialize = 1.1508e-7, doc = 'diffusivity of hydrogen: 1.15e-7 (m^2/s)') # function of T
    m.H_H2 = Param(initialize = 25479, doc = 'Henry constant (bar m^3/mol)')

    # Other parameters
    m.dGstd = Param(initialize = 237.14e3, doc = 'standrad Gibbs free energy(J/mol)')
    m.R = Param(initialize = 0.035e-3, doc = 'Electrode resistance(ohm * m)') # note
    m.rho = Param(initialize = 0.30, doc = 'Resistance(ohm*m)') # use this instead
 
    # Auxilary contants
    m.dz = Param(initialize = m.L/(len(m.n)-1), doc='Discretization step size')  

    m.eps = Param(initialize = 1e-7, doc='Prevent zero division error')
    eps = m.eps
    # Equations
    # ================== Initial conditions ===================

    def _BC1(m, i):
        if i == 'H2O':
            return m.Ca[i,0] == 55.5e3
        return Constraint.Skip
    m.BC1 = Constraint(m.i, rule = _BC1, doc = 'Anode BC')

    def _BC2(m, ci):
        if ci == 'H2O':
            return m.Cc[ci,0] == 20e3
        elif ci == 'H2':
            return m.Cc[ci,0] == 2.16e3
    m.BC2 = Constraint(m.ci, rule = _BC2, doc = 'Cathode BC')

    # ================== Equations ===================

    def _Eq1(m, n):
        """
        MB: H2O anode
        """
        if (n == m.n.first()):
            return Constraint.Skip 
        return (- (m.ua[n]*m.Ca['H2O', n] - m.ua[n-1]*m.Ca['H2O', n-1])/m.dz - \
            1/m.h_ch*(m.Nd[n] + m.Neo[n] + m.Nrxn['H2O', n]) == 0 )
    m.Eq1 = Constraint(m.n, rule=_Eq1, doc= 'MB H2O anode')

    def _Eq2(m, n):
        """
        MB: H2O cathode
        """
        if (n == m.n.first()):
            return Constraint.Skip 
        return - (m.uc[n]*m.Cc['H2O', n] - m.uc[n-1]*m.Cc['H2O', n-1])/m.dz +  \
            1/m.h_ch*(m.Nd[n] + m.Neo[n]) == 0
    m.Eq2 = Constraint(m.n, rule=_Eq2, doc= 'MB H2O cathode')

    def _Eq3(m, n):
        """
        MB: O2 anode
        """
        if (n == m.n.first()):
            return Constraint.Skip    
        return - (m.ua[n]*m.Ca['O2', n] - m.ua[n-1]*m.Ca['O2', n-1])/m.dz + \
            1/m.h_ch*(m.Nrxn['O2', n]) == 0   
    m.Eq3 = Constraint(m.n, rule=_Eq3, doc= 'MB O2 anode')

    def _Eq4(m, n):
        """
        MB: H2 anode
        """
        if (n == m.n.first()):
            return Constraint.Skip 
        return (- (m.ua[n]*m.Ca['H2', n] - m.ua[n-1]*m.Ca['H2', n-1])/m.dz + \
            1/m.h_ch*(m.gamma*m.Nper[n]) == 0 )
    m.Eq4 = Constraint(m.n, rule=_Eq4, doc= 'MB H2 anode')

    def _Eq5(m, n):
        """
        MB: H2 cathode
        """
        if (n == m.n.first()):
            return Constraint.Skip 
        return - (m.uc[n]*m.Cc['H2', n] - m.uc[n-1]*m.Cc['H2', n-1])/m.dz + \
            1/m.h_ch*(m.Nrxn['H2', n] - m.gamma*m.Nper[n]) == 0 
    m.Eq5 = Constraint(m.n, rule = _Eq5, doc = 'MB H2 cathode')

    def _Eq6(m, n):
        """
        MB: defining diffusion of H2O
        """
        return m.Nd[n] == m.D_H2O*(m.Ca['H2O', n] - m.Cc['H2O', n])/m.h
    m.Eq6 = Constraint(m.n, rule = _Eq6, doc = 'Eq diffusion of H2O')

    # Check this later
    def _Eq_ne(m,n):
        """
        Auxilary for electro-osmosis drag coeffficient
        """
        # return m.ne[n] == 0.45
        # Check unit here
        return m.ne[n] == 0.0252*m.Pc[n] - 1.9073*m.j[n] + 0.0189*m.Tm - 2.7892
    m.Eq_ne = Constraint(m.n, rule = _Eq_ne, doc = 'Aux Eq of drag coefficient')

    def _Eq7(m, n):
        """
        Electro-osmosis
        """
        return m.Neo[n] == m.ne[n]*m.j[n]*1e4/m.F 
        # return m.Neo[n] == 0
    m.Eq7 = Constraint(m.n, rule = _Eq7, doc = 'Eq Electro-osmosis')

    def _Eq8(m, n):
        """
        Permeation
        """
        return m.Nper[n] == m.D_H2/m.H_H2*(m.Pc[n]*m.Cc['H2',n]/(m.Cc['H2O',n]+m.Cc['H2',n]+eps) \
            - m.Pa[n]*m.Ca['H2',n]/(m.Ca['H2O',n]+m.Ca['O2',n]+m.Ca['H2',n]+eps))/m.h
    m.Eq8 = Constraint(m.n, rule = _Eq8, doc = 'Eq Permeation of H2')

    def _Eq9(m,n):
        """
        rxn term H2
        """
        return m.Nrxn['H2', n] == m.j[n]*1e4/(2*m.F) - m.gamma*m.Nper[n] # unit conversion 1/cm^2 => 1/m^2
    m.Eq9 = Constraint(m.n, rule = _Eq9, doc = 'Eq Reaction of H2')

    def _Eq10(m,n):
        """
        rxn term H2O
        """
        return m.Nrxn['H2O', n] == m.j[n]*1e4/(2*m.F) - m.gamma*m.Nper[n]  # unit conversion 1/cm^2 => 1/m^2
    m.Eq10 = Constraint(m.n, rule = _Eq10, doc = 'Eq Reaction of H2O')

    def _Eq11(m,n):
        """
        rxn term H2
        """
        return m.Nrxn['O2', n] == m.j[n]*1e4/(4*m.F) - 0.5*m.gamma*m.Nper[n] # unit conversion 1/cm^2 => 1/m^2
    m.Eq11 = Constraint(m.n, rule = _Eq11, doc = 'Eq Reaction of O2')

    def _Eq12(m,n):
        """
        Cell voltage
        """
        return m.V[n] == m.E[n] + m.eta[n]
    m.Eq12 = Constraint(m.n, rule = _Eq12, doc = 'Eq Cell volage balance')

    def _Eq13(m,n):
        """
        Open-circuit voltage 
        """
        # Note that the sign here is confusing, need to ask 
        return m.E[n] == m.dG[n]/(2*m.F)
    m.Eq13 = Constraint(m.n, rule = _Eq13, doc = 'Eq open-circuit voltage')

    def _Eq14(m,n):
        """
        Gibbs free energy
        """
        return  m.dG[n] == m.dGstd
    m.Eq14 = Constraint(m.n, rule = _Eq14, doc = 'Eq Gibbs free energy')        

    def _Eq15(m,n):
        """
        Ohmic potential
        """ 
        return m.eta[n] == m.rho*m.h*m.j[n]*1e4 # 100 for unit m==>cm 
    m.Eq15 = Constraint(m.n, rule = _Eq15, doc = 'Eq Ohmic potential')

    # def _Eq16(m):
    #     """
    #     Specify average current density
    #     """
    #     return sum(m.j[n] for n in m.n) == m.j_avg*(sl+1)
    # m.Eq16 = Constraint(rule = _Eq16, doc = 'average current density')

    m.obj = Objective(expr = 5)

    return m

def init_model(m):
    # Set initial values for concentrations
    m.Ca['H2O', 0].fix(55.5e3)
    m.Cc['H2O', 0].fix(20e3)
    m.Cc['H2', 0].fix(2.16e3)
    m.Ca['H2', 0].fix(0)
    m.Ca['O2', 0].fix(0)
    return 


m = build_memcap()

init_model(m)

opt = SolverFactory('gams')
io_options = dict() 

# io_options['solver'] = "ipopth"
# res = opt.solve(m,
#     tee=True,
#     io_options=io_options)

io_options['solver'] = "baron"
res = opt.solve(m,
    tee=True,
    keepfiles=True,
    symbolic_solver_labels=True,
    #add_options = ['GAMS_MODEL.optfile = 1; option reslim=120; option optcr=0.0;'],
    add_options = ['option reslim=160; option optcr=0.0;'],
    tmpdir='/home/zyuliu/PEMEC/Mv1',
    io_options=io_options)

m.obj.pprint()

# Exporting to excel

l1 = []
l2 = []
l3 = []
col1 = ['H2O', 'O2', 'H2']
df1 = pd.DataFrame(columns = col1)
for v in range(0,21):
    l1.append(value(m.Ca['H2O',v]))
    l2.append(value(m.Ca['H2',v]))
    l3.append(value(m.Ca['O2',v]))
df1['H2O'] = l1
df1['H2'] = l2
df1['O2'] = l3

col2 = ['H2O', 'H2']
df2 = pd.DataFrame(columns = col2)
l1 = []
l2 = []
for v in range(0,21):
    l1.append(value(m.Cc['H2O',v]))
    l2.append(value(m.Cc['H2',v]))
df2['H2O'] = l1
df2['H2'] = l2

col3 = ['H2O', 'H2', 'O2']
df3 = pd.DataFrame(columns = col3)
l1 = []
l2 = []
l3 = []
for v in range(0,21):
    l1.append(value(m.Nrxn['H2O',v]))
    l2.append(value(m.Nrxn['H2',v]))
    l3.append(value(m.Nrxn['O2',v]))
df3['H2O'] = l1
df3['H2'] = l2
df3['O2'] = l3

col4 = ['Nd', 'Neo', 'Nper']
df4 = pd.DataFrame(columns = col4)
l1 = []
l2 = []
l3 = []
for v in range(0,21):
    l1.append(value(m.Nd[v]))
    l2.append(value(m.Neo[v]))
    l3.append(value(m.Nper[v]))
df4['Nd'] = l1
df4['Neo'] = l2
df4['Nper'] = l3

col5 = ['V', 'E', 'ohm']
df5 = pd.DataFrame(columns = col5)
l1 = []
l2 = []
l3 = []
for v in range(0,21):
    l1.append(value(m.V[v]))
    l2.append(value(m.E[v]))
    l3.append(value(m.eta[v]))
df5['V'] = l1
df5['E'] = l2
df5['ohm'] = l3

col6 = ['ne']
df6 = pd.DataFrame(columns = col6)
l1 = []
for v in range(0,21):
    l1.append(value(m.ne[v]))
df6['ne'] = l1

col7 = ['j']
df7 = pd.DataFrame(columns = col7)
l1 = []
for v in range(0,21):
    l1.append(value(m.j[v]))
df7['j'] = l1

with pd.ExcelWriter('Model_V1_output.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='Ca')
    df2.to_excel(writer, sheet_name='Cc')
    df3.to_excel(writer, sheet_name='Nrxn')
    df4.to_excel(writer, sheet_name='fluxes')
    df5.to_excel(writer, sheet_name='potentials')
    df6.to_excel(writer, sheet_name='Auxilary')
    df7.to_excel(writer, sheet_name='Current Density')

print('printing infeasible constraints')

from pyomo.util.infeasible import log_infeasible_constraints
log_infeasible_constraints(m)  

# m.Eq12.pprint()
# m.Neo.pprint()

# for x in range(0,21):
#     print(value(m.ne[x]*m.j[x]*1e4/m.F)) 








        

