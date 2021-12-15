import pandas as pd
from pyomo.environ import *
from datetime import datetime

# ================= Initialization =========================

# Sparser model, better initialization.

def build_model(sl = 30):
    """
    sl = Number of control volumes in the system
    """
    m = ConcreteModel()

    m.sl = Param(initialize = sl)
    m.eps = Param(initialize = 1e-12, doc='Prevent zero division error')
    eps = m.eps

    def V_init(m, n):
        return 2
    def j_init(m, n):
        return 1.2 + (0.98 - 1.2)/sl*n
    def Pa_init(m, n):
        return value(m.Pa0)
    def Pc_init(m):
        return value(m.Pc0)
    def Psat_wt_init(m,n):
        return 0.16     # saturation pressure using Antoine's Eq at 328 K
    def Pa_wt_init(m,n):
        return 0.16     # saturation pressure using Antoine's Eq at 328 K
    def Pa_ox_init(m,n):
        return 3

    m.i = Set(initialize=['H2O', 'O2' ,'H2'], doc='set of species')
    m.ci = Set(within = m.i, initialize = ['H2O', 'H2'], doc='subset of species in cathode')
    # creating sl control volumes
    
    m.n = Set(initialize=[i for i in range(0, sl+1)], ordered=True)

    # ========= Variables =================
    m.Ca = Var(m.i, m.n, domain = NonNegativeReals, doc='Concentration on the anode side (mol/m^3)') # Scaled down 1e3 for H2O
    m.Cc = Var(m.ci, m.n, domain = NonNegativeReals, doc='Concentration on the cathode side(mol/m^3)') # Scaled down 1e3 for H2O and H2
    

    m.Nd = Var(m.n, doc = 'diffusion flux of water(mol/m^2*s)') # Scaled up 1e3
    m.Neo = Var(m.n, doc = 'electro-osmosis flux of water(mol/m^2*s)') # Scaled up 1e2
    m.Nper = Var(m.n, doc = 'Permeation flux of water(mol/m^2*s)') # Scaled up 1e6
    m.Nrxn = Var(m.i, m.n, doc = 'reaction flux(mol/m^2*s)') # Scaled up 1e3 
    
    m.E = Var(m.n, domain = NonNegativeReals, doc = 'open circuit voltage(V)')
    m.dG = Var(m.n, doc = 'Gibbs free energy(J/mol)')
    m.ohm = Var(m.n, domain = NonNegativeReals, doc = 'ohmic resistance(ohm)')
    m.j = Var(m.n, initialize = j_init, domain = NonNegativeReals, doc = 'current density(A/cm^2)')
    m.V = Var(domain = NonNegativeReals, doc = 'cell voltage (V)')
    # m.V = Var(m.n, domain = NonNegativeReals, doc = 'cell voltage (V)')

    # Introduced in v1.1

    m.kappa = Var(m.n, initialize = 0,domain = NonNegativeReals, doc = 'fraction of reaction that is in back permeation')
    m.act = Var(m.n, domain = NonNegativeReals, doc = 'Activation overvoltage')
    m.dif = Var(m.n, domain = NonNegativeReals, doc = 'Diffusion ovorvoltage')

    m.C_O2_wc = Var(m.n, domain = NonNegativeReals, doc = 'O2 anode concentration @ working condition')
    m.C_H2_wc = Var(m.n, domain = NonNegativeReals, doc = 'H2 anode concentration @ working condition')


    # Introduced in 1.6

    def Ca_tot_init(m,n):
        return 55.5e3

    m.Ca_tot = Var(m.n, initialize = Ca_tot_init, doc = 'Total concentration in anode')
    # m.Ca_tot = Var(m.n, doc = 'Total concentration in anode')
    m.Cc_tot = Var(m.n, doc = 'Total concentration in cathode')
    m.cp_a = Var(m.n, doc = 'molar heat capacity in anode')
    m.cp_c = Var(m.n, doc = 'molar heat capacity in cathode')

    # Introduced in 1.8:s

    m.Psat_wt = Var(m.n, initialize = Psat_wt_init, doc = 'Saturation Pressure of Water')
    m.Pa_wt = Var(m.n, initialize = Pa_wt_init, doc = 'Pressure of water in Anode')
    m.Pa_ox = Var(m.n, initialize = Pa_ox_init, doc = 'Pressure of oxygen in the Anode')
    # ========= Auxilary variables =================

    m.ne = Var(m.n, doc = 'Drag coefficient of electro-osmosis(unitless)')
    
    # ========= Parameters ================
    # Constant 
    m.L = Param(initialize=0.3, doc='Lenght of the channel (m)') 
    m.Pa0 = Param(initialize=2, doc='inlet 2 (bar)') # Which pressure to use
    m.Pc0 = Param(initialize=100, doc='inlet 2 (bar)') # Which pressure to use
    m.b = Param(initialize=0.105, doc= '0.105 (m)')
    m.h_ch = Param(initialize = 0.003, doc = 'height of the channel (m)') 
    
    m.h_m = Param(initialize = 0.0002, doc = 'membrane thickness (m)')  
    m.F = Param(initialize = 96485, doc = 'Faradays constant (C/mol)')    
    m.gamma = Param(initialize = 0.2, doc = 'reaction fracction of decomposition(unitless)') 
    m.D_H2O = Param(initialize = 1.28e-10, doc = 'diffusivity of H2O(m^2/s)')
    
    # m.j = Param(m.n, initialize = j_init, doc = 'j_init')
    # m.V = Param(m.n, initialize = V_init, doc = 'v_init')
    # m.V = Param(initialize = 2, doc = 'v_init')
    
    m.Pa = Param(m.n, initialize = Pa_init, doc = 'Anode Pressure (bar) initialization')

    # Coefficient 
    # m.D_H2 = Param(initialize = 1.1508e-7, doc = 'diffusivity of hydrogen: 1.15e-7 (m^2/s)') # function of T
    m.D_H2 = Param(initialize = 1.6297e-6, doc = 'diffusivity of hydrogen: 1.15e-7 (m^2/s)') # function of T
    # m.H_H2 = Param(initialize = 25.479, doc = 'Henry constant (bar m^3/mol)')
    m.H_O2 = Param(initialize = 0.2, doc = 'Henry constant (atm m^3/mol)')
    m.H_H2 = Param(initialize = 25479, doc = 'Henry constant (bar m^3/mol)')
    

    # Other parameters
    m.dGstd = Param(initialize = 237.14e3, doc = 'standrad Gibbs free energy(J/mol)')
    # m.R = Param(initialize = 0.035e-3, doc = 'Electrode resistance(ohm * m)') # note
    m.rho = Param(initialize = 0.30, doc = 'Resistance(ohm*m)') # use this instead
 
    m.R = Param(initialize = 8.3145, doc = 'gas cosntant(J/(mol*K))')
    m.alpha_a = Param(initialize = 2, doc = 'transfer coefficient(anode)') 
    m.alpha_c = Param(initialize = 0.5, doc = 'transfer coefficient(cathode)') 
    m.j0_a = Param(initialize = 10e-10, doc = 'exchange current density at anode(A/cm^2)')
    m.j0_c = Param(initialize = 10e-3, doc = 'exchange current density at cathode(A/cm^2)')
    m.D_H2_bp = Param(initialize = 1.3e-5, doc = 'diffusivity of hydrogen in electode(m^2/s)')
    m.D_O2_bp = Param(initialize = 7.34e-7, doc = 'diffusivity of hydrogen in electode(m^2/s)')
    m.hbp = Param(initialize = 0.0005, doc = 'electode height (m)')  

    # introduced in 1.6

    m.fbp = Param(initialize = 0.1e3, doc='heat transfer coefficient between anode and cathode (W/(K*m2))')
    m.fc = Param(initialize = 0.3e3, doc='heat transfer coefficient between cathode and MEA(W/(K*m2))')
    m.fa = Param(initialize = 0.3e3, doc='heat transfer coefficient between anode and MEA(W/(K*m2))')
    m.cp = Param(m.i, initialize= {'H2O':75.445, 'O2':15.06, 'H2':30.11 }, doc = 'molar heat capacity(J/(mol*K))')
    m.cp_m = Param(initialize = 508, doc = 'specific heat of membrane(J/(kg*K))') 
    m.cp_bp = Param(initialize = 450, doc = 'specific heat of bipolar(J/(kg*K))')
    m.rhom = Param(initialize = 3860, doc='density of MEA(kg/m3)')
    m.rhobp = Param(initialize = 8100, doc = 'density of bipolar plate(kg/m3)')
    m.rhos = Param(m.i, mutable = True, doc = 'density of species')
    m.mw = Param(m.i, mutable = True, doc = 'molecular weight')
    m.Ta_in = Param(initialize = 55 + 273.15, doc = 'Inlet temperature of anode(K)')

    # Auxilary contants
    m.dz = Param(initialize = m.L/(len(m.n)-1), doc='Discretization step size')  
    m.Ca0 = Param(m.i, initialize = {'H2O':55.5e3, 'O2':eps, 'H2':eps }, doc='Initialial conditions')
    m.zH = Param(initialize = 1.05, doc = 'Compressibility')
    m.p2SI = Param(initialize = 1e5, doc = 'converting pressure from bar to Pa') 

   # Initializing parameters

    m.rhos['H2O'] = 1e3 #kg/m^3
    m.mw['H2O'] = 18e-3 #kg/mol

    # More Varibles 
    m.ua_0 = Var(initialize = 2.66, doc = 'velocity(m/min)')  
    m.uc_0 = Param(initialize = 0.144, doc = 'velocity(m/min)') 
    m.uc_l = Param(initialize = 1.2096e-4*60/1.3, doc = 'velocity(m/min)') 
    m.EH = Var(m.n, domain = NonNegativeReals, doc = 'voltage(V) contributing to heat of formation change')

    def ua_init(m, n):
        return value(m.ua_0) + n*(4.8657-m.ua_0)/sl
    def uc_init(m, n):
        # return (value(m.uc_0) + n*(m.uc_l-0.144)/sl)
        return m.uc_l + (n-m.sl)*(m.uc_l-m.uc_0)/sl

    def Ta_init(m,n):
        return 328
    def Tm_init(m,n):
        return 348
    def Tc_init(m,n):
        return 338

    m.Tc_l = Var(initialize = 338, doc = 'Edge case for cathode') 
    m.Ta = Var(m.n, initialize = Ta_init, doc = 'Temperature of the membrane') 
    m.Tm = Var(m.n, initialize = Tm_init, doc = 'Temperature of anode') 
    m.Tc = Var(m.n, initialize = Tc_init, doc = 'Temperature of cathode') 
    # m.Ta = Param(m.n, initialize = Ta_init, doc = 'Temperature of anode') 
    # m.Tm = Param(m.n, initialize = Tm_init, doc = 'Temperature of anode') 
    # m.Tc = Param(m.n, initialize = Tc_init, doc = 'Temperature of cathode') 

    m.ua_ref = Param(m.n, initialize = ua_init, doc = 'Anode velocity' )
    m.uc_ref = Param(m.n, initialize = uc_init, doc = 'Cathode velocity')
    m.ua = Var(m.n, domain = NonNegativeReals, initialize = ua_init, doc = 'Anode velocity' )
    m.uc = Var(m.n, domain = NonNegativeReals, initialize = uc_init, doc = 'Cathode velocity')

    m.j_avg = Var(initialize = 1, doc = 'Average current density')
    # m.j_avg.setlb(0.4)
    # m.j_avg.setub(1.3)
    m.Pc = Var(initialize = Pc_init, doc = 'Cathode Pressure (bar) initialization')
    # m.Pc = Var(initialize = Pc_init, doc = 'Cathode Pressure (bar) initialization')
    # m.Pc.setlb(3)
    # m.Pc.setub(200)

    # ================== Equations ===================

    def _Eq1(m, n):
        """
        MB: H2O anode
        """
        if ( n == m.n.first()):
            return m.Ca['H2O',n] == m.Ca0['H2O'] 
        return ( (m.ua[n]*m.Ca['H2O',n] - m.ua[n-1]*m.Ca['H2O',n-1])/m.dz + \
            60*1/m.h_ch*(m.Nd[n] + m.Neo[n] + m.Nrxn['H2O', n]) == 0 )
    m.Eq1 = Constraint(m.n, rule=_Eq1, doc= 'MB H2O anode')

    def _Eq2(m, n): 
        """
        MB: H2O cathode
        """
        if (n == m.n.last()):
            return ( ( m.uc[n]*m.Cc['H2O',n])/m.dz - \
            60*1/m.h_ch*(m.Nd[n] + m.Neo[n] + m.Nrxn['H2O', n]) == 0 )
        return  ( ( m.uc[n]*m.Cc['H2O',n]-m.uc[n+1]*m.Cc['H2O', n+1])/m.dz - \
            60*1/m.h_ch*(m.Nd[n] + m.Neo[n] + m.Nrxn['H2O', n]) == 0 )
    m.Eq2 = Constraint(m.n, rule=_Eq2, doc= 'MB H2O cathode')

    def _Eq3(m, n):
        """
        MB: O2 Anode
        """
        if ( n == m.n.first()):
            return m.Ca['O2',n] == m.Ca0['O2'] 
        return ( (m.ua[n]*m.Ca['O2', n] - m.ua[n-1]*m.Ca['O2',n-1])/m.dz - \
            60*1/m.h_ch*(m.Nrxn['O2', n]) == 0 )
    m.Eq3 = Constraint(m.n, rule=_Eq3, doc= 'MB O2 anode')

    def _Eq4(m, n):
            """
            MB: H2 Anode
            """
            if ( n == m.n.first()):
                return m.Ca['H2',n] == m.Ca0['H2']
            return ( (m.ua[n]*m.Ca['H2', n] - m.ua[n-1]*m.Ca['H2', n-1])/m.dz - \
                60*1/m.h_ch*(m.gamma*m.Nper[n]) == 0 )
    m.Eq4 = Constraint(m.n, rule=_Eq4, doc= 'MB H2 anode')

    def _Eq5(m, n):
        """
        MB: H2 Cathode
        """
        if (n == m.n.last()):
            return (m.uc[n]*m.Cc['H2', n])/m.dz - \
            60*1/m.h_ch*(m.Nrxn['H2', n] - m.gamma*m.Nper[n]) == 0 
        return (m.uc[n]*m.Cc['H2', n] - m.uc[n+1]*m.Cc['H2', n+1])/m.dz - \
            60*1/m.h_ch*(m.Nrxn['H2', n] - m.gamma*m.Nper[n]) == 0 
    m.Eq5 = Constraint(m.n, rule = _Eq5, doc = 'MB H2 cathode')

    def _Eq6(m, n):
        """
        MB: defining diffusion of H2O
        """
        return m.Nd[n] == m.D_H2O*(m.Ca['H2O', n] - m.Cc['H2O', n])/m.h_m
    m.Eq6 = Constraint(m.n, rule = _Eq6, doc = 'Eq diffusion of H2O')

    def _Eq_ne(m,n):
        """
        Auxilary for electro-osmosis drag coeffficient
        """
        return m.ne[n] == 0.0252*m.Pc - 1.9073*m.j[n] + 0.0189*m.Tm[n] - 2.7892
    m.Eq_ne = Constraint(m.n, rule = _Eq_ne, doc = 'Aux Eq of drag coefficient')

    def _Eq7(m, n):
        """
        Electro-osmosis
        """
        return m.Neo[n] == m.ne[n]*m.j[n]*1e4/m.F 
    m.Eq7 = Constraint(m.n, rule = _Eq7, doc = 'Eq Electro-osmosis')

    def _conc_anode(m, n):
        """
        Sum of concentration of anode
        """
        return m.Ca_tot[n] == m.Ca['H2O',n] + m.Ca['H2',n] + m.Ca['O2',n]
    m.Eq_conc_anode = Constraint(m.n, rule = _conc_anode, doc = 'total concentration of anode')

    def _conc_cathode(m, n):
        """
        Sum of concentration of cathode
        """
        return m.Cc_tot[n] == m.Cc['H2O',n] + m.Cc['H2',n]
    m.Eq_conc_cathode = Constraint(m.n, rule = _conc_cathode, doc = 'total concentration of cathode')
 
    def _Eq8(m, n):
        """
        Permeation
        """
        # return m.Nper[n] == m.D_H2/m.H_H2*(m.Pc-m.Pa[n]*m.Ca['H2',n]/(m.Ca['H2',n]+m.Ca['O2',n]+eps))/m.h_m
        return m.Nper[n] == m.D_H2/m.H_H2*(m.Pc)/m.h_m
    m.Eq8 = Constraint(m.n, rule = _Eq8, doc = 'Eq Permeation of H2')

    def _Eq9(m,n):
        """
        rxn term H2
        """
        return m.Nrxn['H2', n] == m.j[n]*1e4/(2*m.F) - (1-m.gamma)*m.Nper[n] # unit conversion 1/cm^2 => 1/m^2
    m.Eq9 = Constraint(m.n, rule = _Eq9, doc = 'Eq Reaction of H2')

    def _Eq10(m,n):
        """
        rxn term H2O
        """
        return m.Nrxn['H2O', n] == m.j[n]*1e4/(2*m.F) - (1-m.gamma)*m.Nper[n]  # unit conversion 1/cm^2 => 1/m^2
    m.Eq10 = Constraint(m.n, rule = _Eq10, doc = 'Eq Reaction of H2O')

    def _Eq11(m,n):
        """
        rxn term H2
        """
        return m.Nrxn['O2', n] == m.j[n]*1e4/(4*m.F) - 0.5*(1-m.gamma)*m.Nper[n] # unit conversion 1/cm^2 => 1/m^2
    m.Eq11 = Constraint(m.n, rule = _Eq11, doc = 'Eq Reaction of O2')

    def _Eq12(m,n):
        """
        Cell voltage
        """
        return m.V == m.E[n] + m.ohm[n] + m.act[n] + m.dif[n]
        # return m.V == m.E[n] + m.ohm[n] + m.dif[n]
    m.Eq12 = Constraint(m.n, rule = _Eq12, doc = 'Eq Cell volage balance')

    def _Eq13(m,n):
        """
        Open-circuit voltage 
        """
        # Note that the sign here is confusing, need to ask 
        return m.E[n] == (1-m.kappa[n])*m.dG[n]/(2*m.F)
    m.Eq13 = Constraint(m.n, rule = _Eq13, doc = 'Eq open-circuit voltage')

    def _Eq14(m,n):
        """
        Kappa
        """
        return m.kappa[n] == (1-m.gamma)*m.Nper[n]/(m.Nrxn['H2',n]+((1-m.gamma)*m.Nper[n])+eps)
    m.Eq14 = Constraint(m.n, rule = _Eq14, doc = 'Eq fraction of back permeation')

    def _Eqwt_Antoine(m,n): 
        return m.Psat_wt[n] == 0.158
    m.Eqwt_Antoine = Constraint(m.n, rule = _Eqwt_Antoine, doc = 'Saturation Pressure of Water in Anode')

    def _Eqwt_Pa(m,n):
        return m.Pa_wt[n] ==  (m.Ca['H2O',n]/(m.Ca_tot[n]+eps))*m.Psat_wt[n]
    m.Eqwt_Pa = Constraint(m.n, rule = _Eqwt_Pa, doc = 'Partial Pressure of Water in Anode')

    def _Eqox_Pa(m,n):
        return m.Pa_ox[n] == m.H_O2*m.Ca['O2',n]
    m.Eqox_Pa = Constraint(m.n, rule = _Eqox_Pa, doc = 'Pressure at Eq of O2 in Anode')

    def _Eq15(m,n):
        """
        Gibbs free energy
        """
        # return  m.dG[n] == m.dGstd + m.R*m.Tm[n]*log(m.Pc) - m.R*m.Tm[n]*log(m.Pa_wt[n]+eps)\
        #     + m.R*m.Tm[n]*log(eps+(-1.864*2.71828**(-0.518*m.Ca['O2',n]) + 2.0385)**0.5)
        return  m.dG[n] == m.dGstd + m.R*m.Tm[n]*log(m.Pc) + 0.5*m.R*m.Tm[n]*log(0.2+m.Pa[n]*m.Pa_ox[n]/(m.Pa_ox[n]+m.Pa_wt[n])) - \
            m.R*m.Tm[n]*log(eps+m.Pa[n]*m.Pa_wt[n]/(m.Pa_ox[n]+m.Pa_wt[n]))
    m.Eq15 = Constraint(m.n, rule = _Eq15, doc = 'Eq Gibbs free energy') 

    def _Eq16(m,n):
        """
        Ohmic potential
        """ 
        return m.ohm[n] == m.rho*m.h_m*m.j[n]*1e4 # 100 for unit m==>cm 
    m.Eq16 = Constraint(m.n, rule = _Eq16, doc = 'Eq Ohmic potential')

    def arcsinh(x):
        return log(x+(x**2+1)**0.5)

    def arcsinh(x):
        return log(x+(x**2+1)**0.5)

    def _Eq17(m,n):
        """
        Activation potential
        """
        # return m.act[n] == m.R*m.Ta[n]/(m.alpha_a*m.F)*arcsinh(m.j[n]/(2*m.j0_a)) +\
        #     m.R*m.Tc[n]/(m.alpha_c*m.F)*arcsinh(m.j[n]/(2*m.j0_c))
        return m.act[n] == m.R*m.Ta[n]/(m.alpha_a*m.F)*1.44 +\
            m.R*m.Tc[n]/(m.alpha_c*m.F)*0.001
    m.Eq17 = Constraint(m.n, rule = _Eq17, doc = 'Eq activation potential')

    def _Eq18(m,n):
        """
        Difussion potential
        """
        return m.dif[n] == m.R*m.Ta[n]/(4*m.F) + m.R*m.Tc[n]/(2*m.F)
        # return m.dif[n] == m.R*m.Ta[n]/(4*m.F)*log(eps + m.Ca['O2',n]/eps + m.C_O2_wc[n]) + m.R*m.Tc[n]/(2*m.F)*log(eps + m.Cc['H2',n]/eps + m.C_H2_wc[n])
    m.Eq18 = Constraint(m.n, rule = _Eq18, doc = 'Eq diffusion potential')

    def _Eq19(m,n):
        """
        working condition concentration of H2
        """
        return m.Nrxn['H2',n] == m.D_H2_bp*(m.Cc['H2',n]-m.C_H2_wc[n])/m.hbp
    m.Eq19 = Constraint(m.n, rule = _Eq19, doc = 'Eq C_H2 wc')

    def _Eq20(m,n):
        """
        working condition concentration of O2
        """
        return m.Nrxn['O2',n] == m.D_O2_bp*(m.C_O2_wc[n]-m.Ca['O2',n])/m.hbp
    m.Eq20 = Constraint(m.n, rule = _Eq20, doc = 'Eq C_O2 wc')

    def _Eq21(m,n):
        """
        Molar heat capacity in anode
        """
        return m.cp_a[n] == (m.Ca['H2O',n]*m.cp['H2O'] + m.Ca['O2',n]*m.cp['O2']\
         + m.Ca['H2',n]*m.cp['H2'])/(m.Ca_tot[n]+eps)
    m.Eq21 = Constraint(m.n, rule = _Eq21, doc = 'Molar heat capacity(Anode)')

    def _Eq22(m,n):
        """
        Molar heat capacity in cathode
        """
        return m.cp_c[n] == (m.Cc['H2O',n]*m.cp['H2O'] + m.Cc['H2',n]*m.cp['H2'])/(m.Cc_tot[n]+eps)
    m.Eq22 = Constraint(m.n, rule = _Eq22, doc = 'Molar heat capacity(Cathode)')

    def _Eq_EH(m,n):
        """
        Thermo-neutral voltage
        """
        # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.852.9751&rep=rep1&type=pdf
        return m.EH[n] == 1.48 
    m.Eq_EH = Constraint(m.n, rule = _Eq_EH, doc = 'Thermoneutral voltage')

    def _Eq_jreq(m):
        return sum(m.j[n] for n in m.n) == m.j_avg*(sl+1)
    m.Eq_jreq = Constraint(rule = _Eq_jreq, doc = 'Controlling average current density')

    m.obj = Objective(expr = 5)

    return m

def init_model(m):
    m.Cc['H2O', m.n.last()] = 19.91e3
    m.Cc['H2', m.n.last()] = 2.14e3
    m.hhv = Param(initialize = 285.8e3, doc = 'HHV for the combustion of hydrogen(J/mol)')
    
    m.Qf = Var(initialize = 5/120, within = NonNegativeReals, doc = 'L/min')
    m.Q0 = Var(initialize = 100/120, within = NonNegativeReals, doc = 'L/min')
    m.Qi = Var(initialize = 95/100, within = NonNegativeReals)
    m.Pc.fix(value(m.Pc))
    m.j_avg.fix(value(m.j_avg))
    # m.ua_0.fix(value(m.ua_0))
    m.Q0.fix(value(m.Q0))

    def _Eq_Qi(m):
        return m.Qi + m.Qf == m.Q0
    m.Eq_Qi = Constraint(rule = _Eq_Qi, doc = 'feed flow rate') 

    def _Eq_Q0(m):
        return (m.Q0)*1e-3/(m.b*m.h_ch) == m.ua_0 
    m.Eq_Q0 = Constraint(rule = _Eq_Q0, doc = 'recycle balance')

    def _Eq_Qf(m):
        return m.Qf == m.uc[0]*(m.b*m.h_ch)*1e3 
    m.Eq_Qf = Constraint(rule = _Eq_Qf, doc = 'outlet flowrate')

    for np in m.n:
        m.ua[np].fix(m.ua_ref[np])
        m.uc[np].fix(m.uc_ref[np])
        m.Ta[np].fix(328)
        m.Tc[np].fix(338)
        m.Tm[np].fix(348)
    m.Tc_l.fix(338)

    # m.j.setlb(0.7)
    return 

def add_energy_balance(m):
    for np in m.n:
        m.Ta[np].unfix()
        m.Tc[np].unfix()
        m.Tm[np].unfix()
        m.Tc_l.unfix()
    eps = m.eps

    def _Eq23(m,n):
        """
        EB: Anode
        """
        if (n == m.n.first()):
            return m.Ta[n] == m.Ta_in
        return m.ua[n]*(m.Ta[n]-m.Ta[n-1])/m.dz == 60*(m.fbp/(m.Ca_tot[n]*m.cp_a[n]*m.h_ch+eps)*(m.Tc[n]-m.Ta[n])+\
            (m.cp['H2']*m.gamma*m.Nper[n] + m.cp['O2']*m.Nrxn['O2',n]+m.fa)*(m.Tm[n]-m.Ta[n])/(m.Ca_tot[n]*\
            m.cp_a[n]*m.h_ch+eps))
    m.Eq23 = Constraint(m.n, rule = _Eq23, doc = 'Energy Balance(Anode)')

    def _Eq_Tcl(m):
        return m.Tc_l == m.Ta[value(m.sl)] + m.Tc[0]-m.Ta[0]
    m.Eq_Tcl = Constraint(rule = _Eq_Tcl)


    def _Eq24(m,n):
        """
        EB: Cathode
        """
        if (n == m.n.last()):
            return m.uc[n]*(m.Tc[n]-m.Tc_l)/m.dz == 60*(m.fbp/(m.Cc_tot[n]*m.cp_c[n]*m.h_ch+eps)*(m.Ta[n]-m.Tc[n])+\
                (m.cp['H2O']*(m.Nd[n] + m.Neo[n]) + m.cp['H2']*(m.Nrxn['H2',n]+(1-m.gamma)*m.Nper[n]) +m.fc)*(m.Tm[n]-m.Tc[n])/(m.Cc_tot[n]*\
                m.cp_c[n]*m.h_ch+eps))
        return m.uc[n]*(m.Tc[n]-m.Tc[n+1])/m.dz == 60*(m.fbp/(m.Cc_tot[n]*m.cp_c[n]*m.h_ch+eps)*(m.Ta[n]-m.Tc[n])+\
            (m.cp['H2O']*(m.Nd[n] + m.Neo[n]) + m.cp['H2']*(m.Nrxn['H2',n]+(1-m.gamma)*m.Nper[n]) +m.fc)*(m.Tm[n]-m.Tc[n])/(m.Cc_tot[n]*\
            m.cp_c[n]*m.h_ch+eps))
    m.Eq24 = Constraint(m.n, rule = _Eq24, doc = 'Energy Balance(Cathode)')

    def _Eq25(m,n):
        """
        EB: MEA-bipolar
        """
        return (m.V-m.EH[n])*m.j[n]*1e4/(m.h_m*m.rhom*m.cp_m + m.hbp*m.rhobp*m.cp_bp) + (m.cp['H2O']*(m.Nd[n]+\
            m.Neo[n]+m.Nrxn['H2O',n])+m.fa)*(m.Ta[n]-m.Tm[n])/(m.h_m*m.rhom*m.cp_m + m.hbp*m.rhobp*m.cp_bp) +\
            (m.cp['H2']*m.Nper[n]+m.fc)*(m.Tc[n]-m.Tm[n])/(m.h_m*m.rhom*m.cp_m + m.hbp*m.rhobp*m.cp_bp) == 0
    m.Eq25 = Constraint(m.n, rule = _Eq25, doc = 'Energy Balance(MEA-bipolar)')

    return 

def expand_model(m):
    for np in m.n:
        m.ua[np].unfix()
        m.uc[np].unfix()
 
    def _Eq_ua_dat(m,n):
        if n == m.n.first():
            return m.ua[0] == m.ua_0
        return m.ua_0 + n*m.j_avg*0.0367617*60/m.sl == m.ua[n]
    m.Eq_ua_dat = Constraint(m.n, rule = _Eq_ua_dat)

    def _Eq_uc_dat(m,n):
       return m.uc_l + (n-m.sl)*(m.uc_l-m.uc_0)/m.sl == m.uc[n]
    m.Eq_uc_dat = Constraint(m.n, rule = _Eq_uc_dat)

    return 

def set_velocity(m):
    m.del_component(m.Eq_ua_dat)
    m.del_component(m.Eq_uc_dat)

    def _Eq_ua(m,n):
        if n == m.n.first():
            return m.ua[0] == m.ua_0
        return 1 == m.Ca['H2O',n]*m.mw['H2O']/m.rhos['H2O'] + (m.Ca['O2',n]+m.Ca['H2',n])*m.R*m.Ta[n]/(m.Pa[n]*100000)
    m.Eq_ua = Constraint(m.n, rule = _Eq_ua)

    def _Eq_uc(m,n):
        return 1 == m.Cc['H2O',n]*m.mw['H2O']/(m.rhos['H2O']*0.976) + (m.Cc['H2',n])*m.R*m.Tc[n]/(m.zH*m.Pc*100000)
    m.Eq_uc = Constraint(m.n, rule = _Eq_uc)

    return 

def set_current_req(m, javg=1):
    m.del_component(m.Eq_jreq)
    m.j_avg.fix(javg)
    def _Eq_jreq(m):
        return sum(m.j[n] for n in m.n) == m.j_avg*(m.sl+1)
    m.Eq_jreq = Constraint(rule = _Eq_jreq, doc = 'Controlling average current density')
    return 

def add_sweep(m, objf=0, relax=0):
    """
        prepare for parametric sweep for different obj
    """
    # m.del_component(m.obj)

    # m.ecost = Param(initialize = 6.66e-2, doc = '$/kWh')
    m.ecost = Param(initialize = 13.89e-2, doc = '$/kWh') # pitt avg
    m.Hprice = Param(initialize = 3, doc = '$/kg' )

    m.Hprod = Var(doc = 'Hydrogen Production(mol/s)')
    m.P_op = Var(doc = 'power used in operating condition')
    m.eta_sys = Var(doc = 'system efficiency')
    m.cost = Var(doc = 'cost at operating condition')
    m.Hprod_price = Var(doc = 'Hydrogen produce in dollar value')

    def _Eq_Hprod(m):
        """ mol/s """
        return m.Hprod == m.b*m.h_ch*m.Cc['H2',0]*m.uc[0]
        # return m.Hprod == m.hhv*m.Cc['H2',0]*m.uc[0]
    m.Eq_Hprod = Constraint(rule =  _Eq_Hprod)

    def _Eq_P_op(m):
        """ W """
        return m.P_op == 1e4*m.j_avg*m.V*m.L*m.b
        # return m.P_op == 1e4*m.j_avg*m.V
    m.Eq_P_op = Constraint(rule = _Eq_P_op)

    def _Eq_eta_sys(m):
        """ % """
        return  m.eta_sys == m.hhv*m.Hprod/(m.P_op+m.eps)*100
        #  return 1e-3*m.eta_sys == (1-m.hhv*m.Hprod/(m.P_op+m.eps))*(m.P_op+m.eps)/m.hhv*18e-3*m.Hprice*3600
    m.Eq_eta_sys = Constraint(rule = _Eq_eta_sys)  

    def _Eq_cost(m):
        """ 1e-3 $/h """
        return m.cost == m.ecost*m.P_op
    m.Eq_cost = Constraint(rule = _Eq_cost)

    def _Eq_Hprod_price(m):
        """ 1e-3$/h """
        return 1e-3*m.Hprod_price == m.Hprod*18e-3*3600*m.Hprice
    m.Eq_Hprod_price = Constraint(rule = _Eq_Hprod_price)

    return

m = build_model()


init_model(m)

opt = SolverFactory('gams')
io_options = dict() 

now = datetime.now()
start_time = now.strftime("%H:%M:%S")

io_options['solver'] = "baron"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=60; option optcr=0.0;'],
    io_options=io_options)

io_options['solver'] = "ipopt"
res = opt.solve(m,
    tee=True,
    # keepfiles=True,
    add_options = ['option reslim=90; option optcr=0.0;'],
    # tmpdir='/home/zyuliu/PEMEC/Mv1',
    io_options=io_options)

expand_model(m)

io_options['solver'] = "baron"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=60; option optcr=0.0;'],
    io_options=io_options)

io_options['solver'] = "ipopt"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=90; option optcr=0.0;'],
    io_options=io_options)

add_energy_balance(m)

io_options['solver'] = "baron"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=60; option optcr=0.0;'],
    io_options=io_options)

io_options['solver'] = "ipopt"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=90; option optcr=0.0;'],
    io_options=io_options)

set_velocity(m)

io_options['solver'] = "baron"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=60; option optcr=0.0;'],
    io_options=io_options)

io_options['solver'] = "ipopt"
res = opt.solve(m,
    tee=True,
    add_options = ['option reslim=60; option optcr=0.0;'],
    io_options=io_options)

# Paramter Sweep
# j_list = [0.2, 0.6, 0.8, 1, 1.2, 1.4]
# j_list = [1]
# rep_list = []
# for javg in j_list:
#     set_current_req(m, javg)

#     io_options['solver'] = "ipopt"
#     res = opt.solve(m,
#         tee=True,
#         # keepfiles=True,
#         add_options = ['option reslim=90; option optcr=0.0;'],
#         # tmpdir='/home/zyuliu/PEMEC/Mv1',
#         io_options=io_options)
#     rep_list.append(value(-100*m.hhv*m.Cc['H2',0]*m.uc[0]/(1e4*m.j_avg*m.V)))

# set_current_req(m, 1.2)

io_options['solver'] = "ipopt"
res = opt.solve(m,
    tee=True,
    keepfiles=True,
    # symbolic_solver_labels=True,
    add_options = ['GAMS_MODEL.optfile = 1; option reslim=90; option optcr=0.0;'],
    tmpdir='/home/zyuliu/PEMEC/Mv1',
    io_options=io_options)

now = datetime.now()
end_time = now.strftime("%H:%M:%S")

# io_options['solver'] = "conopt"
# res = opt.solve(m,
#     tee=True,
#     add_options = ['option reslim=90; option optcr=0.0;'],
#     io_options=io_options)

# m.obj.pprint()

# Exporting to excel

sl = 30

l1 = []
l2 = []
l3 = []
col1 = ['H2O', 'O2', 'H2']
df1 = pd.DataFrame(columns = col1)
for v in range(0,sl+1):
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
for v in range(0,sl+1):
    l1.append(value(m.Cc['H2O',v]))
    l2.append(value(m.Cc['H2',v]))
df2['H2O'] = l1
df2['H2'] = l2

col3 = ['H2O', 'H2', 'O2']
df3 = pd.DataFrame(columns = col3)
l1 = []
l2 = []
l3 = []
for v in range(0,sl+1):
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
for v in range(0,sl+1):
    l1.append(value(m.Nd[v]))
    l2.append(value(m.Neo[v]))
    l3.append(value(m.Nper[v]))
df4['Nd'] = l1
df4['Neo'] = l2
df4['Nper'] = l3

col5 = ['V', 'E', 'ohm', 'activation', 'diffusion']
df5 = pd.DataFrame(columns = col5)
l1 = []
l2 = []
l3 = []
l4 = []
l5 = []
for v in range(0,sl+1):
    # l1.append(value(m.V[v]))
    l1.append(value(m.V))
    l2.append(value(m.E[v]))
    l3.append(value(m.ohm[v]))
    l4.append(value(m.act[v]))
    l5.append(value(m.dif[v]))
df5['V'] = l1
df5['E'] = l2
df5['ohm'] = l3
df5['activation'] = l4
df5['diffusion'] = l5

col6 = ['ne']
df6 = pd.DataFrame(columns = col6)
l1 = []
for v in range(0,sl+1):
    l1.append(value(m.ne[v]))
df6['ne'] = l1

col7 = ['j']
df7 = pd.DataFrame(columns = col7)
l1 = []
for v in range(0,sl+1):
    l1.append(value(m.j[v]))
df7['j'] = l1

col8 = ['Ca_O2_wc','Cc_H2_wc']
df8 = pd.DataFrame(columns = col8)
l1 = []
l2 = []
for v in range(0,sl+1):
    l1.append(value(m.C_O2_wc[v]))
    l2.append(value(m.C_H2_wc[v]))
df8['Ca_O2_wc'] = l1 
df8['Cc_H2_wc'] = l2

col9 = ['Ua','Uc']
df9 = pd.DataFrame(columns = col9)
l1 = []
l2 = []
l1.append(value(m.ua_0)*60)
l2.append(value(m.uc_0)*60)
for v in range(0,sl+1):
    l1.append(value(m.ua[v])*60)
    l2.append(value(m.uc[v])*60)
df9['Ua'] = l1
df9['Uc'] = l2

col10 = ['Ta','Tc','Tm']
df10 = pd.DataFrame(columns = col10)
l1 = []
l2 = []
l3 = []
for v in range(0,sl+1):
    l1.append(value(m.Ta[v])-273.15)
    l2.append(value(m.Tc[v])-273.15)
    l3.append(value(m.Tm[v])-273.15)
df10['Ta'] = l1
df10['Tc'] = l2
df10['Tm'] = l3

with pd.ExcelWriter('Model_V2_0_output.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='Ca')
    df2.to_excel(writer, sheet_name='Cc')
    df3.to_excel(writer, sheet_name='Nrxn')
    df4.to_excel(writer, sheet_name='fluxes')
    df5.to_excel(writer, sheet_name='potentials')
    df6.to_excel(writer, sheet_name='Auxilary')
    df7.to_excel(writer, sheet_name='Current Density')
    df8.to_excel(writer, sheet_name='Working conditions')
    df9.to_excel(writer, sheet_name='Velocity')
    df10.to_excel(writer, sheet_name='Temperature')

# print('H2 production')
# for i in rep_list:
#     print(i)
# print(value(1e3*m.Cc['H2',0]*m.uc[0]))

print('printing infeasible constraints')

from pyomo.util.infeasible import log_infeasible_constraints
log_infeasible_constraints(m)


print(value(m.Qf)*120)
print(value(m.Q0)*120)
print(value(m.kappa[15]))
print(start_time)
print(end_time)
