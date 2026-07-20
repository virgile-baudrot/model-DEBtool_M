function [par, metaPar, txtPar] = pars_init_Turdus_philomelos_CT36EFSA(metaData)

metaPar.model = 'std'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 9000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 2.3308;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.04527;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.98119;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 510.7866;   free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.019945;   free.k_J   = 1;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 7364.1283;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 1.709e+02; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hp = 3.758e+03; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.220e-12;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.001;      free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hx = 1.153e+03; free.E_Hx  = 1;   units.E_Hx = 'J';         label.E_Hx = 'maturity at fledging'; 
par.f_embryo = 1;       free.f_embryo = 0;    units.f_embryo = '-';   label.f_embryo = 'scaled functional response during embryo stage';

% --- Normal Baselines (Memory for starvation windows) ---
par.f_base_nest = 2.0;        free.f_base_nest = 0;  units.f_base_nest = '-';  label.f_base_nest = 'baseline f during nest';
par.f_base_juvenile  = 1.0;   free.f_base_juvenile  = 0;  units.f_base_juvenile  = '-';  label.f_base_juvenile  = 'baseline f during juvenile';
par.f_base_adult   = 1.0;     free.f_base_adult   = 0;  units.f_base_adult   = '-';  label.f_base_adult   = 'baseline f during adult';

% --- Active parameters (These get modified by the runner script) ---
par.f_nest     = par.f_base_nest;       free.f_nest     = 0;  units.f_nest     = '-';  label.f_nest     = 'active f during nest';
par.f_juvenile = par.f_base_juvenile;   free.f_juvenile = 0;  units.f_juvenile = '-';  label.f_juvenile = 'active f during juvenile';
par.f_adult    = par.f_base_adult;      free.f_adult    = 0;  units.f_adult    = '-';  label.f_adult    = 'active f during adult';

par.t_R = 365;        free.t_R   = 0;   units.t_R = 'd';          label.t_R = 'time since birth at first brood'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free;