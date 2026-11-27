function [prdData, info] = predict_Podarcis_muralis_CT36EFSA(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  % compute temperature correction factors
  TC_tp = tempcorr(temp.tp, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_ab = tempcorr(C2K(Tab(:,1)), T_ref, T_A);
  TC_K  = tempcorr(temp.tL_K, T_ref, T_A);
  TC_D  = tempcorr(temp.tL_D, T_ref, T_A);
  TC_LN = tempcorr(temp.LN, T_ref, T_A);
  
  % Male rescaling properties
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap * E_m_m);            % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
  L_mm = v/ k_M/ g_m;                  % cm, max struct length

  % ---------- zero-variate data: CHRONIC outputs
  % Each stage f is held constant from birth (tp at f_juvenile) or over the
  % whole life (Li, Wwi, Ri, am at f_adult); f_embryo is the mother's reserve
  % density at egg production (egg size, reserve at hatching).
  % An event that cannot be reached at that f (ultimate maturity below the
  % threshold, see get_tp.m) is returned as NaN; the run itself goes on.
  % (The male g_m only rescales the parameters: same condition with v_Hpm.)
  reach = @(f_val, vH) f_val * max(0, f_val - l_T)^2 > k * vH; % -, DEBtool condition

  % life cycle
  pars_tp = [g k l_T v_Hb v_Hp];
  [t_b, l_b, info] = get_tb([g; k; v_Hb], f_embryo); % -, scaled age & length at birth
  if ~info; prdData = []; return; end
  [t_p, t_bp, l_p] = event_CT36EFSA(pars_tp, f_juvenile, reach(f_juvenile, v_Hp));
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; 
  U_E0 = initial_scaled_reserve(f_embryo, pars_UE0); 
  Ww_0 = U_E0 * p_Am * w_E/ mu_E/ d_E; 
  
  % birth
  L_b = L_m * l_b;                  
  Lw_b = L_b/ del_M;                
  Ww_b = L_b^3 * (1 + f_embryo * w);       

  % puberty 
  L_p = L_m * l_p;                  
  Lw_p = L_p/ del_M;                
  tT_p = (t_p-t_bp)/ k_M/ TC_tp;    % d, time since birth at puberty (NaN if not reached)

  % ultimate
  l_i = f_adult - l_T;                    
  L_i = L_m * l_i;                  
  Lw_i = L_i/ del_M;                
  Ww_i = L_i^3 * (1 + f_adult * w);       

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; 
  if reach(f_adult, v_Hp)
    RT_i = TC_Ri * reprod_rate(L_i, f_adult, pars_R);
  else
    RT_i = 0; % maturity maintenance at E_Hp exceeds the flux to maturity/reproduction
  end

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  
  t_m = get_tm_s(pars_tm, f_adult, l_b);      
  aT_m = t_m/ k_M/ TC_am;               
  
  % males events
  L_im = L_mm * f_adult;
  Lw_im = L_im/ del_Mm; 
  pars_tpm = [g_m k l_T v_Hb v_Hpm];
  l_bm = l_b * z/ z_m;                 % -, male scaled length at birth (same egg as females)
  [~, ~, l_pm] = event_CT36EFSA(pars_tpm, f_juvenile, reach(f_juvenile, v_Hpm), [t_b, f_juvenile, l_bm]);
  Lw_pm = L_mm * l_pm/ del_Mm;
  Ww_im = L_im^3 * (1 + f_adult * w_m);      

  % pack to output
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Lpm = Lw_pm;
  prdData.Li = Lw_i;
  prdData.Lim = Lw_im;
  prdData.Ww0 = Ww_0;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Wwim = Ww_im;
  prdData.Ri = RT_i;
  
  % ---------- uni-variate data: ACUTE scenario (n-day window)
  %% --- CONTINUOUS ODE INTEGRATION (Robust against starvation freezes) ---
  d=7; % days, duration of altered f values
  % Stage onsets = BASELINE events (unrestricted f), so that the window starts
  % at the onset of the stage and the baseline timeline is the same for every
  % reduction factor.
  [t_p0, t_bp0] = get_tp(pars_tp, f_base_juvenile);
  [t_pm0, t_bpm0] = get_tp(pars_tpm, f_base_juvenile, [t_b, f_base_juvenile, l_bm]);
  tT_p = (t_p0 - t_bp0)/ k_M/ TC_tp;    % d, baseline female time since birth at puberty
  % Full-period scenario (APPENDIX C, par.ct36_full = 1): the reduced f is kept over the
  % whole stage, which then ends at the event reached at the reduced f (Inf if never reached).
  is_full = isfield(par, 'ct36_full') && par.ct36_full;
  if is_full
    d = Inf;
    if f_juvenile ~= f_base_juvenile; tT_p = (t_p - t_bp)/ k_M/ TC_tp; tT_p(isnan(tT_p)) = Inf; end
  end
  % Time vector provided from runner script
  t_run = t_run; 
  tau_run_K = t_run * k_M * TC_K; % For plots matching 20C data
  tau_run_D = t_run * k_M * TC_D; 

  % === FEMALE ===
  Tf = unique([0; tT_p; min(d, tT_p); tT_p + d]); Tf = Tf(Tf >= 0 & isfinite(Tf)); % d, exact switch times of f
  f_vec_f = f_base_juvenile * ones(size(Tf));
  f_vec_f(Tf >= tT_p) = f_base_adult;
  
  if f_juvenile ~= f_base_juvenile
      f_vec_f(Tf >= 0 & Tf < min(d, tT_p)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_f(Tf >= tT_p & Tf < tT_p + d) = f_adult;
  end
  
  
  el_0_f = [1; l_b]; 
  el_res_fK = el_path_CT36EFSA(tau_run_K, Tf * k_M * TC_K, f_vec_f, el_0_f, g, l_T); % exact n-day windows
  el_res_fD = el_path_CT36EFSA(tau_run_D, Tf * k_M * TC_D, f_vec_f, el_0_f, g, l_T); % exact n-day windows
  
  ELw_fK = L_m .* el_res_fK(:, 2) / del_M;
  EWw_fK = (L_m .* el_res_fK(:, 2)).^3 .* (1 + el_res_fK(:, 1) .* w);
  ELw_fD = L_m .* el_res_fD(:, 2) / del_M;

  % === MALE ===
  tT_pm = (t_pm0 - t_bpm0)/ k_M/ TC_tp; % d, baseline male time since birth at puberty
  if is_full && f_juvenile ~= f_base_juvenile % male puberty at the reduced f
    [t_pm1, t_bpm1] = event_CT36EFSA(pars_tpm, f_juvenile, reach(f_juvenile, v_Hpm), [t_b, f_juvenile, l_bm]);
    tT_pm = (t_pm1 - t_bpm1)/ k_M/ TC_tp; tT_pm(isnan(tT_pm)) = Inf;
  end
  Tm = unique([0; tT_pm; min(d, tT_pm); tT_pm + d]); Tm = Tm(Tm >= 0 & isfinite(Tm)); % d, exact switch times of f
  f_vec_m = f_base_juvenile * ones(size(Tm));
  f_vec_m(Tm >= tT_pm) = f_base_adult;

  
  if f_juvenile ~= f_base_juvenile
      f_vec_m(Tm >= 0 & Tm < min(d, tT_pm)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_m(Tm >= tT_pm & Tm < tT_pm + d) = f_adult;
  end
  
  
  el_0_m = [1; l_bm]; 
  el_res_mK = el_path_CT36EFSA(tau_run_K, Tm * k_M * TC_K, f_vec_m, el_0_m, g_m, l_T); % exact n-day windows
  el_res_mD = el_path_CT36EFSA(tau_run_D, Tm * k_M * TC_D, f_vec_m, el_0_m, g_m, l_T); % exact n-day windows
  
  ELw_mK = L_mm .* el_res_mK(:, 2) / del_Mm;
  EWw_mK = (L_mm .* el_res_mK(:, 2)).^3 .* (1 + el_res_mK(:, 1) .* w_m);
  ELw_mD = L_mm .* el_res_mD(:, 2) / del_Mm;

  %% Extract exact points required by standard DEBtool 'mydata'
  prdData.tL_K = [interp1(t_run, ELw_fK, tL_K(:,1), 'linear', 'extrap'), ...
                  interp1(t_run, ELw_mK, tL_K(:,1), 'linear', 'extrap')];
  prdData.tL_D = [interp1(t_run, ELw_fD, tL_D(:,1), 'linear', 'extrap'), ...
                  interp1(t_run, ELw_mD, tL_D(:,1), 'linear', 'extrap')];
              
  prdData.LW_f = interp1(ELw_fK, EWw_fK, LW_f(:,1), 'linear', 'extrap');
  prdData.LW_m = interp1(ELw_mK, EWw_mK, LW_m(:,1), 'linear', 'extrap');
  
  prdData.LN = TC_LN * 365 * reprod_rate(LN(:,1) * del_M, f_adult, pars_R);
  prdData.Tab = t_b/ k_M ./ TC_ab;
  
  %% Provide continuous points for beautiful plotting in runner script
  prdData.plot_t = t_run;
  prdData.plot_SVL_f = ELw_fK;
  prdData.plot_SVL_m = ELw_mK;
  prdData.plot_Ww_f  = EWw_fK;
  prdData.plot_Ww_m  = EWw_mK;
  prdData.plot_L_run = L_run;
  if reach(f_adult, v_Hp)
    prdData.plot_LN    = TC_LN * 365 * reprod_rate(L_run * del_M, f_adult, pars_R);
  else
    prdData.plot_LN = 0 * L_run; % no reproduction below the chronic threshold (as Ri)
  end

end

% Local custom ODE function for robust reserve (e) and length (l) integration
function del = dget_el_local(tau, el, tf_mat, g, l_T)
    e = max(1e-6, el(1)); 
    l = max(1e-6, el(2)); 

    f = interp1(tf_mat(:,1), tf_mat(:,2), tau, 'previous', 'extrap');

    de = g * (f - e) / l;
    rho = g * (e/l - 1 - l_T/l) / (e + g);
    dl = l * rho / 3;

    del = [de; dl];
end