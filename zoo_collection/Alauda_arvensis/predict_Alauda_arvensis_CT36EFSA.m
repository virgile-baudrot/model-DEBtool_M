function [prdData, info] = predict_Alauda_arvensis_CT36EFSA(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;

  % zero-variate data
  % life cycle
  pars_tx = [g; k; l_T; v_Hb; v_Hx];  % compose parameter vector
  [tau_x, ~, l_x, ~, info] = get_tp(pars_tx, f_nest);
  pars_tp = [g; k; l_T; v_Hb; v_Hp];  % compose parameter vector
  [~, tau_b, ~, l_b, info] = get_tp(pars_tp, f_embryo);  % -, scaled times & lengths
  [tau_p, ~, l_p, ~, info] = get_tp(pars_tp, f_juvenile);  % -, scaled times & lengths
  if ~info; prdData = []; return; end

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * ome);       % g, wet weight at birth at f 
  tT_b = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T
  aT_b = tT_b;                % d, time at birth at f and T
  t_g = tT_b;                 % d, gestation time (time from fertilization to birth)
  
  % fledging
  tT_x = (tau_x - tau_b)/ kT_M;         % d, time since birth at fledging
  L_x = L_m * l_x;
  Ww_x = L_x^3 * (1 + f_nest * ome);

  % puberty
  tT_p = (tau_p - tau_b)/ kT_M;         % d, time since birth at puberty
  L_p = L_m * l_p;
  Ww_p = L_p^3 * (1 + f_juvenile * ome);

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * ome);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate(L_i, f_adult, pars_R);                % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  tau_m = get_tm_s(pars_tm, f_adult, l_b);      % -, scaled mean life span at T_ref
  aT_m = tau_m/ kT_M;                     % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tg = t_g;
  prdData.tx = tT_x;
  prdData.tp = tT_p;
  prdData.tR = t_R;
  prdData.am = aT_m;
  prdData.Li = L_i;
  prdData.Wwb = Ww_b;
  prdData.Wwx = Ww_x;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  
  % uni-variate data
  % time-weight: varying food
  % tf = [3 6 9 12; f_03 f_06 f_09 f_12]'; tf(:,1) = tf(:,1) * kT_M;
  % tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1)*kT_M);
  % EWw = (L_m * tvel(:,4)).^3 .* (1 + tvel(:,3) * ome);   % g, weight

  %% --- NEW PIECEWISE FOOD WITH A 7-DAY RESTRICTED STARVATION WINDOW ---
  % 1. Build the unstarved baseline timeline for the whole life cycle
  % (Using the baseline parameters defined in pars_init)
  f = f_base_nest * ones(size(tW(:,1)));
  f(tW(:,1) >= tT_x) = f_base_juvenile;
  f(tW(:,1) >= tT_p) = f_base_adult;
  % 2. Apply the 7-day starvation ONLY if the runner modified the parameter
  % Nest stage starvation (starts at birth t=0, lasts 7 days)
  if f_nest ~= f_base_nest
      f(tW(:,1) >= 0 & tW(:,1) < min(7, tT_x)) = f_nest;
  end
  % Juvenile stage starvation (starts at fledging t=tT_x, lasts 7 days)
  if f_juvenile ~= f_base_juvenile
      f(tW(:,1) >= tT_x & tW(:,1) < tT_x + 7) = f_juvenile;
  end
  % Adult stage starvation (starts at puberty t=tT_p, lasts 7 days)
  if f_adult ~= f_base_adult
      f(tW(:,1) >= tT_p & tW(:,1) < tT_p + 7) = f_adult;
  end
  % Build time–food table expected by get_tp
  tf = [tW(:,1) * kT_M, f];

  % Call get_tp once (vectorized in time, piecewise in f)
  tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1) * kT_M);
  % Compute wet weight
  EWw = (L_m .* tvel(:,4)).^3 .* (1 + tvel(:,3) .* ome);

  % pack to output
  prdData.tW = EWw;

end