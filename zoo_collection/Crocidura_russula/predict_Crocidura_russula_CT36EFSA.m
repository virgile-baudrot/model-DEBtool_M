function [prdData, info] = predict_Crocidura_russula_CT36EFSA(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  % compute temperature correction factors
  TC = tempcorr(temp.tg, T_ref, T_A);
  
  % zero-variate data

  % life cycle
  pars_tx = [g k l_T v_Hb v_Hx v_Hp]; info = 1;
  [~, ~, t_b, ~, ~, l_b] = get_tx(pars_tx, f_embryo);
  [t_p, t_x, ~, l_p, l_x, ~] = get_tx(pars_tx, f_juvenile);
 
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * w);       % g, wet weight at birth at f 
  tT_g = t_0 + t_b/ k_M/ TC;        % d, gestation time at f and T

  % weaning
  L_x = L_m * l_x;                  % cm, structural length at weaning
  Ww_x = L_x^3 * (1 + f_nest * w);       % g, wet weight at weaning
  tT_x = (t_x - t_b)/ k_M/ TC;      % d, time since birth at fledging at f and T
  
  % puberty 
  tT_p = (t_p - t_b)/ k_M/ TC;      % d, time since birth at puberty at f and T

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate_foetus(L_i, f_adult, pars_R);         % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_adult, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;                  % d, mean life span at T
  
  % pack to output
  prdData.tg = tT_g;
  prdData.tx = tT_x;
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Wwb = Ww_b;
  prdData.Wwx = Ww_x;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  
  % uni-variate data
  
  % time-weight 
  f = f_nest * ones(size(tW(:,1)));   % initialize with nest (under weaning) values
  f(tW(:,1) >= tT_x) = f_juvenile;        % replace where time >= tT_x
  f(tW(:,1) >= tT_p) = f_adult;           % replace where time >= tT_p

  %{
  pars_lb = [g; k; v_Hb];                 % parameters
  ir_B = 3/k_M + 3 .* f .* L_m ./ v;      % d, von Bert growth rate
  rT_B = TC ./ ir_B; 
    
  L_b = get_lb(pars_lb, f_embryo) * L_m;     % cm, length at birth
  L_i = f .* L_m - L_T;                      % cm, ultimate length
  L   = L_i - (L_i - L_b) .* exp(- rT_B .* tW(:,1)); % cm, structural length at time
  EWw = L.^3 .* (1 + f .* w);                % g, wet weight
  %}
  % Build time–food table expected by get_tp
  kT_M=1;
  pars_lb = [g; k; v_Hb];
  pars_tb = [g; k; v_Hb];
  % scaled length at birth
  [l_b, info_lb] = get_lb(pars_lb, f_embryo); if info_lb == 0 prdData = []; return; end
  [tau_b, info_tb] = get_tb(pars_tb, f_embryo); if info_tb == 0 prdData = []; return; end

  pars_tp = [g; k; l_T; v_Hb; v_Hp];  % compose parameter vector
  tf = [tW(:,1) * kT_M, f];
  % Call get_tp once (vectorized in time, piecewise in f)
  tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1) * kT_M);
  % Compute wet weight
  EWw = (L_m .* tvel(:,4)).^3 .* (1 + tvel(:,3) .* ome);
  
  % pack to output
  prdData.tW = EWw;
  