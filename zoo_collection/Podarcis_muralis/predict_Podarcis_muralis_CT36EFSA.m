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

  %% --- PRE-CHECK: Prevent ODE45 freezing for both sexes ---
  v_H_max   = @(f_val) f_val * max(0, f_val - l_T)^2 * (max(0, f_val - l_T) + g) / (f_val + g) / k;
  v_H_max_m = @(f_val) f_val * max(0, f_val - l_T)^2 * (max(0, f_val - l_T) + g_m) / (f_val + g_m) / k;
  
  if v_H_max(f_embryo) < v_Hp || v_H_max(f_juvenile) < v_Hp || v_H_max(f_adult) < v_Hp
      info = 0; prdData = []; return;
  end
  if v_H_max_m(f_embryo) < v_Hpm || v_H_max_m(f_juvenile) < v_Hpm || v_H_max_m(f_adult) < v_Hpm
      info = 0; prdData = []; return;
  end
  %% --------------------------------------------------------

  % life cycle
  pars_tp = [g k l_T v_Hb v_Hp];
  [~, t_b, ~, l_b, info] = get_tp(pars_tp, f_embryo);     
  if ~info; prdData = []; return; end
  [t_p, ~, l_p, ~, ~] = get_tp(pars_tp, f_juvenile);
  
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
  tT_p = (t_p-t_b)/ k_M/ TC_tp;     

  % ultimate
  l_i = f_adult - l_T;                    
  L_i = L_m * l_i;                  
  Lw_i = L_i/ del_M;                
  Ww_i = L_i^3 * (1 + f_adult * w);       

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; 
  RT_i = TC_Ri * reprod_rate(L_i, f_adult, pars_R);             

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  
  t_m = get_tm_s(pars_tm, f_adult, l_b);      
  aT_m = t_m/ k_M/ TC_am;               
  
  % males events
  L_im = L_mm * f_adult;
  Lw_im = L_im/ del_Mm; 
  pars_tpm = [g_m k l_T v_Hb v_Hpm];
  [~, t_bm, ~, l_bm] = get_tp(pars_tpm, f_embryo, [l_b f_embryo*z/z_m l_b*z/z_m]);
  [t_pm, ~, l_pm, ~] = get_tp(pars_tpm, f_juvenile, [l_b f_juvenile*z/z_m l_b*z/z_m]);
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
  
  %% --- CONTINUOUS ODE INTEGRATION (Robust against starvation freezes) ---
  d=7; % days, duration of altered f values
  % Time vector provided from runner script
  t_run = t_run; 
  tau_run_K = t_run * k_M * TC_K; % For plots matching 20C data
  tau_run_D = t_run * k_M * TC_D; 

  % === FEMALE ===
  f_vec_f = f_base_juvenile * ones(size(t_run));
  f_vec_f(t_run >= tT_p) = f_base_adult;
  
  if f_juvenile ~= f_base_juvenile
      f_vec_f(t_run >= 0 & t_run < min(d, tT_p)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_f(t_run >= tT_p & t_run < tT_p + d) = f_adult;
  end
  
  tf_mat_fK = [[-1; tau_run_K; tau_run_K(end)+100], [f_vec_f(1); f_vec_f; f_vec_f(end)]];
  tf_mat_fD = [[-1; tau_run_D; tau_run_D(end)+100], [f_vec_f(1); f_vec_f; f_vec_f(end)]];
  
  el_0_f = [1; l_b]; 
  [~, el_res_fK] = ode45(@(t, y) dget_el_local(t, y, tf_mat_fK, g, l_T), tau_run_K, el_0_f);
  [~, el_res_fD] = ode45(@(t, y) dget_el_local(t, y, tf_mat_fD, g, l_T), tau_run_D, el_0_f);
  
  ELw_fK = L_m .* el_res_fK(:, 2) / del_M;
  EWw_fK = (L_m .* el_res_fK(:, 2)).^3 .* (1 + el_res_fK(:, 1) .* w);
  ELw_fD = L_m .* el_res_fD(:, 2) / del_M;

  % === MALE ===
  tT_pm = (t_pm - t_bm)/ k_M/ TC_tp;
  f_vec_m = f_base_juvenile * ones(size(t_run));
  f_vec_m(t_run >= tT_pm) = f_base_adult;

  
  if f_juvenile ~= f_base_juvenile
      f_vec_m(t_run >= 0 & t_run < min(d, tT_pm)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec_m(t_run >= tT_pm & t_run < tT_pm + d) = f_adult;
  end
  
  tf_mat_mK = [[-1; tau_run_K; tau_run_K(end)+100], [f_vec_m(1); f_vec_m; f_vec_m(end)]];
  tf_mat_mD = [[-1; tau_run_D; tau_run_D(end)+100], [f_vec_m(1); f_vec_m; f_vec_m(end)]];
  
  el_0_m = [1; l_bm]; 
  [~, el_res_mK] = ode45(@(t, y) dget_el_local(t, y, tf_mat_mK, g_m, l_T), tau_run_K, el_0_m);
  [~, el_res_mD] = ode45(@(t, y) dget_el_local(t, y, tf_mat_mD, g_m, l_T), tau_run_D, el_0_m);
  
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
  prdData.plot_LN    = TC_LN * 365 * reprod_rate(L_run * del_M, f_adult, pars_R);

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