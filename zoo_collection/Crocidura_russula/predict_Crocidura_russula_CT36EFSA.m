function [prdData, info] = predict_Crocidura_russula_CT36EFSA(par, data, auxData)

  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  % compute temperature correction factors
  TC = tempcorr(temp.tg, T_ref, T_A);
  % zero-variate data
  pars_tx = [g k l_T v_Hb v_Hx v_Hp];
  %% --- PRE-CHECK: Prevent ODE45 freezing ---
  v_H_max = @(f_val) f_val * max(0, f_val - l_T)^2 * (max(0, f_val - l_T) + g) / (f_val + g) / k;
  % AJOUT DE f_adult DANS LE FILTRE DE SÉCURITÉ
  if v_H_max(f_foetus) < v_Hp || v_H_max(f_lactation) < v_Hp || v_H_max(f_juvenile) < v_Hp || v_H_max(f_adult) < v_Hp
      info = 0; prdData = []; return;
  end
  %% -----------------------------------------

  % life cycle - Piecewise computation for stx mammal model
  [~, ~, tau_b, ~, ~, l_b, info] = get_tx(pars_tx, f_foetus);
  [~, tau_x, ~, ~, l_x, ~] = get_tx(pars_tx, f_lactation);
  [tau_p, ~, ~, l_p, ~, ~] = get_tx(pars_tx, f_juvenile);
  if ~info; prdData = []; return; end

  % birth
  L_b = L_m * l_b;                           % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_foetus * w);         % g, wet weight at birth at f 
  tT_b = tau_b / k_M / TC;                   % d, time from fertilization to birth
  tT_g = t_0 + tT_b;                         % d, gestation time at f and T

  % weaning
  L_x = L_m * l_x;                           % cm, structural length at weaning
  Ww_x = L_x^3 * (1 + f_lactation * w);      % g, wet weight at weaning
  tT_x = (tau_x - tau_b)/ k_M / TC;          % d, time since birth at weaning at f and T

  % puberty
  L_p = L_m * l_p;                           % cm, structural length at puberty at f
  Ww_p = L_p^3 *(1 + f_juvenile * w);        % g, wet weight at puberty 
  tT_p = (tau_p - tau_b)/ k_M / TC;          % d, time since birth at puberty at f and T

  % ultimate
  l_i = f_adult - l_T;                       % -, scaled ultimate length
  L_i = L_m * l_i;                           % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * w);          % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate_foetus(L_i, f_adult, pars_R);   % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];       % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_adult, l_b);     % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;                       % d, mean life span at T
  % pack to output
  prdData.tg = tT_g;
  prdData.tx = tT_x;
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Wwb = Ww_b;
  prdData.Wwx = Ww_x;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  % uni-variate data (Time-Weight)

  %% --- NEW PIECEWISE FOOD WITH A 7-DAY RESTRICTED STARVATION WINDOW ---
  % 1. Build the unstarved baseline timeline for the whole life cycle
  f_vec = f_base_lactation * ones(size(tW(:,1)));
  f_vec(tW(:,1) >= tT_x) = f_base_juvenile;
  f_vec(tW(:,1) >= tT_p) = f_base_adult;

  % 2. Apply the 7-day starvation ONLY if the runner modified the parameter
  if f_lactation ~= f_base_lactation
      f_vec(tW(:,1) >= 0 & tW(:,1) < min(7, tT_x)) = f_lactation;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec(tW(:,1) >= tT_x & tW(:,1) < tT_x + 7) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec(tW(:,1) >= tT_p & tW(:,1) < tT_p + 7) = f_adult;
  end

  % Time axis for post-natal development in DEBtool is scaled time since birth
  tau_run = tW(:,1) * k_M * TC;
  tf_mat = [[-1; tau_run; tau_run(end)+100], [f_vec(1); f_vec; f_vec(end)]];

  % Custom robust ODE integration bypassing get_tx event detection
  el_0 = [1; l_b]; % Reserve density (e) at birth is 1, scaled length is l_b
  [~, el_res] = ode45(@(t, y) dget_el_local(t, y, tf_mat, g, l_T), tau_run, el_0);

  % Calculate wet weight from the simulated structures and reserves
  L_t = L_m .* el_res(:, 2);
  EWw = L_t.^3 .* (1 + el_res(:, 1) .* w);

  % pack to output
  prdData.tW = EWw;

end

% Local custom ODE function for reserve (e) and length (l)
function del = dget_el_local(tau, el, tf_mat, g, l_T)
    e = max(1e-6, el(1)); 
    l = max(1e-6, el(2)); 

    f = interp1(tf_mat(:,1), tf_mat(:,2), tau, 'previous', 'extrap');

    de = g * (f - e) / l;
    rho = g * (e/l - 1 - l_T/l) / (e + g);
    dl = l * rho / 3;
    del = [de; dl];
end