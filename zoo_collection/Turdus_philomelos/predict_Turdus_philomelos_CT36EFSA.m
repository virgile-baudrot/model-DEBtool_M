function [prdData, info] = predict_Turdus_philomelos(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;

  % zero-variate data

  % life cycle
  pars_tx = [g; k; l_T; v_Hb; v_Hx];  % compose parameter vector
  tau_x = get_tp(pars_tx, f_nest); 
  pars_tp = [g; k; l_T; v_Hb; v_Hp];  % compose parameter vector
  [~, tau_b, ~, l_b, info] = get_tp(pars_tp, f_embryo);  % -, scaled times & lengths
  [tau_p, ~, l_p, ~, info] = get_tp(pars_tp, f_juvenile);  % -, scaled times & lengths
  if ~info; prdData = []; return; end

  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * ome);       % g, wet weight at birth at f 
  aT_b = tau_b/ k_M/ TC_ab;     % d, age at birth at f and T

  % fledging
  tT_x = (tau_x - tau_b)/ kT_M;         % d, time since birth at fledging

  % puberty
  tT_p = (tau_p - tau_b)/ kT_M;         % d, time since birth at puberty

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * ome);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate(L_i, f_adult, pars_R);             % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  tau_m = get_tm_s(pars_tm, f_adult, l_b);      % -, scaled mean life span at T_ref
  aT_m = tau_m/ kT_M;                     % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tx = tT_x;
  prdData.tp = tT_p;
  prdData.tR = t_R;
  prdData.am = aT_m;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  
  % uni-variate data
  % time-weight
  % tf = [0 03 05 07 09 10 12; f_00 f_03 f_05 f_07 f_09 f_10 f_12]'; tf(:,1) = tf(:,1) * kT_M;
  % tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1)*kT_M);
  % EWw = (L_m * tvel(:,4)).^3 .* (1 + tvel(:,3) * ome);   % g, weight

  % Build piecewise-constant f(t)
  f = f_nest * ones(size(tW(:,1)));
  f(tW(:,1) >= tT_x) = f_juvenile;
  f(tW(:,1) >= tT_p) = f_adult;
  % Build time–food table expected by get_tp
  tf = [tW(:,1) * kT_M, f];
  % Call get_tp once (vectorized in time, piecewise in f)
  tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1) * kT_M);
  % Compute wet weight
  EWw = (L_m .* tvel(:,4)).^3 .* (1 + tvel(:,3) .* ome);
  
  % pack to output
  prdData.tW = EWw;

end