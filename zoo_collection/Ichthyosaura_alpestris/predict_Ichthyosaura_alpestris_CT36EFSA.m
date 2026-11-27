function [prdData, info] = predict_Ichthyosaura_alpestris_CT36EFSA(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  if E_Hpm < E_Hb
    prdData = []; info = 0; return
  end

  % compute temperature correction factors
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;
  
  % Male rescaling properties
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
  L_mm = v/ k_M/ g_m;                  % cm, max struct length

  % ---------- zero-variate data: CHRONIC outputs
  % Each stage f is held constant from birth (tj at f_nest, tp at f_juvenile)
  % or over the whole life (Li, Wwi, Ri, am at f_adult).
  % An event that cannot be reached at that f (ultimate maturity below the
  % threshold, see get_tp.m) is returned as NaN; the run itself goes on.
  % (The male g_m only rescales the parameters: same condition with v_Hpm.)
  reach = @(f_val, vH) f_val * max(0, f_val - l_T)^2 > k * vH; % -, DEBtool condition

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj];               % compose parameter vector
  pars_tp = [g; k; l_T; v_Hb; v_Hp];               % compose parameter vector
  [tau_b, l_b, info] = get_tb([g; k; v_Hb], f_embryo); % -, scaled age & length at birth
  if ~info; prdData = []; return; end
  [tau_j, tau_bj, l_j] = event_CT36EFSA(pars_tj, f_nest, reach(f_nest, v_Hj));
  [tau_p, tau_bp, l_p] = event_CT36EFSA(pars_tp, f_juvenile, reach(f_juvenile, v_Hp));
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * w);       % g, wet weight at birth at f
  
  % metam
  L_j = L_m * l_j;                  % cm, struc length at metam
  Lw_j = L_j/ del_M;                % cm, SVL at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, struc length at puberty
  Lw_p = L_p/ del_M;                % cm, SVL at puberty

  % metam
  tT_j = (tau_j - tau_bj)/ kT_M;         % d, time since birth at metamorphosis (NaN if not reached)
  % puberty
  tT_p = (tau_p - tau_bp)/ kT_M;         % d, time since birth at puberty (NaN if not reached)

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate SVL at f
  Ww_i = L_i^3 * (1 + f_adult * w);       % g, ultimate wet weight

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  if reach(f_adult, v_Hp)
    RT_i = TC * reprod_rate(L_i, f_adult, pars_R);% #/d, reprod rate
  else
    RT_i = 0; % maturity maintenance at E_Hp exceeds the flux to maturity/reproduction
  end

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_adult, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ kT_M;                     % d, mean life span at T
  
  % Male events
  pars_tpm = [g_m; k; l_T; v_Hb; v_Hpm]; pars_tjm = [g_m; k; l_T; v_Hb; v_Hj];
  [tau_bm, l_bm] = get_tb([g_m; k; v_Hb], f_embryo);
  [tau_pm, ~, l_pm] = event_CT36EFSA(pars_tpm, f_juvenile, reach(f_juvenile, v_Hpm));

  L_pm = L_mm * l_pm; Lw_pm = L_pm/ del_M;   % cm, SVL at puberty
  L_im = f_adult * L_mm; Lw_im = L_im/ del_M;% cm, ultimate struct, SVL
  Ww_im = L_im^3 * (1 + f_adult * w_m);      % g, ultimate wet weight

  % pack to output
  prdData.am = aT_m;
  prdData.tj = tT_j;
  prdData.tp = tT_p;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Lpm = Lw_pm;
  prdData.Li = Lw_i;
  prdData.Lim = Lw_im;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Wwim = Ww_im;
  prdData.Ri = RT_i;
  
  % ---------- uni-variate data: ACUTE scenario (n-day window)
  %% --- PIECEWISE FOOD WITH A 336-DAY RESTRICTED WINDOW ---
  d_nest= 7*4*12;
  d = 7*4*12; % d, 48 weeks = 336 days
  % Stage onsets = BASELINE events (unrestricted f), so that the window starts
  % at the onset of the stage and the baseline timeline is the same for every
  % reduction factor.
  [tau_j0, tau_bj0] = get_tp(pars_tj, f_base_nest);
  [tau_p0, tau_bp0] = get_tp(pars_tp, f_base_juvenile);
  tT_j = (tau_j0 - tau_bj0)/ kT_M;       % d, baseline time since birth at metamorphosis
  tT_p = (tau_p0 - tau_bp0)/ kT_M;       % d, baseline time since birth at puberty
  [tau_jm0, tau_bjm0] = get_tp(pars_tjm, f_base_nest);
  [tau_pm0, tau_bpm0] = get_tp(pars_tpm, f_base_juvenile);
  tT_jm = (tau_jm0 - tau_bjm0)/ kT_M;    % d, baseline male time since birth at metamorphosis
  tT_pm = (tau_pm0 - tau_bpm0)/ kT_M;    % d, baseline male time since birth at puberty
  % Full-period scenario (APPENDIX C, par.ct36_full = 1): the reduced f is kept over the
  % whole stage, which then ends at the event reached at the reduced f (Inf if never reached).
  is_full = isfield(par, 'ct36_full') && par.ct36_full;
  if is_full
    d = Inf; d_nest = Inf;
    if f_nest ~= f_base_nest; tT_j = (tau_j - tau_bj)/ kT_M; tT_j(isnan(tT_j)) = Inf; end
    if f_juvenile ~= f_base_juvenile; tT_p = (tau_p - tau_bp)/ kT_M; tT_p(isnan(tT_p)) = Inf; end
  end
  if is_full % male events at the reduced f
    if f_nest ~= f_base_nest
      [tau_jm1, tau_bjm1] = event_CT36EFSA(pars_tjm, f_nest, reach(f_nest, v_Hj));
      tT_jm = (tau_jm1 - tau_bjm1)/ kT_M; tT_jm(isnan(tT_jm)) = Inf;
    end
    if f_juvenile ~= f_base_juvenile
      [tau_pm1, tau_bpm1] = event_CT36EFSA(pars_tpm, f_juvenile, reach(f_juvenile, v_Hpm));
      tT_pm = (tau_pm1 - tau_bpm1)/ kT_M; tT_pm(isnan(tT_pm)) = Inf;
    end
  end
  % === FEMALE ===
  Tf = unique([0; tT_j; tT_p; min(d_nest, tT_j); tT_j + d; tT_p + d]); Tf = Tf(Tf >= 0 & isfinite(Tf)); % d, exact switch times of f
  f_vec_f = f_base_nest * ones(size(Tf));
  f_vec_f(Tf >= tT_j) = f_base_juvenile;
  f_vec_f(Tf >= tT_p) = f_base_adult;
  
  if f_nest ~= f_base_nest
      f_vec_f(Tf >= 0 & Tf < min(d_nest, tT_j)) = f_nest;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec_f(Tf >= tT_j & Tf < min(tT_j + d, tT_p)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_f(Tf >= tT_p & Tf < tT_p + d) = f_adult;
  end
  
  tau_run_f = tL_f(:,1) * kT_M;
  
  el_0_f = [1; l_b]; 
  el_res_f = el_path_CT36EFSA(tau_run_f, Tf * kT_M, f_vec_f, el_0_f, g, l_T); % exact n-day windows
  ELw_f = L_m .* el_res_f(:, 2) / del_M;
  
  % === MALE ===
  Tm = unique([0; tT_jm; tT_pm; min(d_nest, tT_jm); tT_jm + d; tT_pm + d]); Tm = Tm(Tm >= 0 & isfinite(Tm)); % d, exact switch times of f
  f_vec_m = f_base_nest * ones(size(Tm));
  f_vec_m(Tm >= tT_jm) = f_base_juvenile;
  f_vec_m(Tm >= tT_pm) = f_base_adult;
  
  if f_nest ~= f_base_nest
      f_vec_m(Tm >= 0 & Tm < min(d_nest, tT_jm)) = f_nest;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec_m(Tm >= tT_jm & Tm < min(tT_jm + d, tT_pm)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_m(Tm >= tT_pm & Tm < tT_pm + d) = f_adult;
  end
  
  tau_run_m = tL_m(:,1) * kT_M;
  
  el_0_m = [1; l_bm]; 
  el_res_m = el_path_CT36EFSA(tau_run_m, Tm * kT_M, f_vec_m, el_0_m, g_m, l_T); % exact n-day windows
  ELw_m = L_mm .* el_res_m(:, 2) / del_M;

  % pack to output
  prdData.tL_f = ELw_f;
  prdData.tL_m = ELw_m;
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