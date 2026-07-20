function [prdData, info] = predict_Rana_temporaria_CT36EFSA(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;
  
  % Male rescaling properties
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
  L_mm = v/ k_M/ g_m;                  % cm, max struct length

  %% --- PRE-CHECK: Prevent ODE45 freezing for both sexes ---
  v_H_max   = @(f_val) f_val * max(0, f_val - l_T)^2 * (max(0, f_val - l_T) + g) / (f_val + g) / k;
  v_H_max_m = @(f_val) f_val * max(0, f_val - l_T)^2 * (max(0, f_val - l_T) + g_m) / (f_val + g_m) / k;
  
  if v_H_max(f_embryo) < v_Hp || v_H_max(f_nest) < v_Hp || v_H_max(f_juvenile) < v_Hp || v_H_max(f_adult) < v_Hp
      info = 0; prdData = []; return;
  end
  if v_H_max_m(f_embryo) < v_Hpm || v_H_max_m(f_nest) < v_Hpm || v_H_max_m(f_juvenile) < v_Hpm || v_H_max_m(f_adult) < v_Hpm
      info = 0; prdData = []; return;
  end
  %% --------------------------------------------------------

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj];               % compose parameter vector
  [~, tau_b, ~, l_b] = get_tp(pars_tj, f_embryo);   % -, scaled times & lengths at f
  [tau_j, ~, l_j, ~] = get_tp(pars_tj, f_nest);     % -, scaled times & lengths at f
  pars_tp = [g; k; l_T; v_Hb; v_Hp];               % compose parameter vector
  [tau_p, ~, l_p, ~, info] = get_tp(pars_tp, f_juvenile); 
  if ~info; prdData = []; return; end
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * w);       % g, wet weight at birth at f
  aT_b = tau_b / k_M / TC_ab;
  
  % metam
  L_j = L_m * l_j;                  % cm, struc length at metam
  Lw_j = L_j/ del_M;                % cm, SVL at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, struc length at puberty
  Lw_p = L_p/ del_M;                % cm, SVL at puberty

  % metam
  tT_j = (tau_j - tau_b)/ kT_M;          % d, time since birth at metamorphosis
  % puberty 
  tT_p = (tau_p - tau_b)/ kT_M;          % d, time since birth at puberty

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate SVL at f
  Ww_i = L_i^3 * (1 + f_adult * w);       % g, ultimate wet weight

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate(L_i, f_adult, pars_R);% #/d, reprod rate 

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_adult, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ kT_M;                     % d, mean life span at T
  
  % Male events
  pars_tpm = [g_m; k; l_T; v_Hb; v_Hpm]; pars_tjm = [g_m; k; l_T; v_Hb; v_Hj];
  [~, tau_bm, ~, l_bm] = get_tp(pars_tjm, f_embryo);
  [tau_jm, ~, ~, ~] = get_tp(pars_tjm, f_nest);
  [tau_pm, ~, l_pm, ~] = get_tp(pars_tpm, f_juvenile);
  
  L_pm = L_mm * l_pm; Lw_pm = L_pm/ del_M;   % cm, SVL at puberty
  L_im = f_adult * L_mm; Lw_im = L_im/ del_M;% cm, ultimate struct, SVL
  Ww_im = L_im^3 * (1 + f_adult * w_m);      % g, ultimate wet weight

  % pack to output
  prdData.ab = aT_b;
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
  
  %% --- NEW PIECEWISE FOOD WITH A 7-DAY RESTRICTED STARVATION WINDOW ---
  
  % === FEMALE TIME LENGTH ===
  f_vec_f = f_base_nest * ones(size(tL_f(:,1)));
  f_vec_f(tL_f(:,1) >= tT_j) = f_base_juvenile;
  f_vec_f(tL_f(:,1) >= tT_p) = f_base_adult;

  d = 7*3*6; % 7 days * weeks
  
  if f_nest ~= f_base_nest
      f_vec_f(tL_f(:,1) >= 0 & tL_f(:,1) < min(d, tT_j)) = f_nest;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec_f(tL_f(:,1) >= tT_j & tL_f(:,1) < min(tT_j + d, tT_p)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_f(tL_f(:,1) >= tT_p & tL_f(:,1) < tT_p + d) = f_adult;
  end
  
  tau_run_f = tL_f(:,1) * kT_M;
  tf_mat_f = [[-1; tau_run_f; tau_run_f(end)+100], [f_vec_f(1); f_vec_f; f_vec_f(end)]];
  
  el_0_f = [1; l_b]; 
  [~, el_res_f] = ode45(@(t, y) dget_el_local(t, y, tf_mat_f, g, l_T), tau_run_f, el_0_f);
  ELw_f = L_m .* el_res_f(:, 2) / del_M;
  
  % === MALE TIME LENGTH ===
  tT_jm = (tau_jm - tau_bm)/ kT_M;
  tT_pm = (tau_pm - tau_bm)/ kT_M;

  f_vec_m = f_base_nest * ones(size(tL_m(:,1)));
  f_vec_m(tL_m(:,1) >= tT_jm) = f_base_juvenile;
  f_vec_m(tL_m(:,1) >= tT_pm) = f_base_adult;
  
  if f_nest ~= f_base_nest
      f_vec_m(tL_m(:,1) >= 0 & tL_m(:,1) < min(d, tT_jm)) = f_nest;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec_m(tL_m(:,1) >= tT_jm & tL_m(:,1) < min(tT_jm + d, tT_pm)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_m(tL_m(:,1) >= tT_pm & tL_m(:,1) < tT_pm + d) = f_adult;
  end
  
  tau_run_m = tL_m(:,1) * kT_M;
  tf_mat_m = [[-1; tau_run_m; tau_run_m(end)+100], [f_vec_m(1); f_vec_m; f_vec_m(end)]];
  
  el_0_m = [1; l_bm]; 
  [~, el_res_m] = ode45(@(t, y) dget_el_local(t, y, tf_mat_m, g_m, l_T), tau_run_m, el_0_m);
  ELw_m = L_mm .* el_res_m(:, 2) / del_M;

  % pack to output
  prdData.tL_f = ELw_f;
  prdData.tL_m = ELw_m;

  % === UNIVARIATE LENGTH-CHANGE (LdL) ===
  % We use the specific f_LdL parameter to match the environmental conditions 
  % in which these specific animals were measured, preserving the realistic slope.
  rT_B = TC_ab * k_M/ 3/ (1 + f_LdL*f_adult/ g); 
  L_i_dL = f_LdL*f_adult * L_m; 
  EdLw = rT_B * (L_i_dL/del_M - LdL(:,1));
  prdData.LdL = EdLw;

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