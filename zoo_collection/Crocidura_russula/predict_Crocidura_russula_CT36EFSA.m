function [prdData, info] = predict_Crocidura_russula_CT36EFSA(par, data, auxData)

  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  % compute temperature correction factors
  TC = tempcorr(temp.tg, T_ref, T_A);
  % ---------- zero-variate data: CHRONIC outputs
  % Each stage f is held constant from birth (tx at f_lactation, tp at f_juvenile)
  % or over the whole life (Li, Wwi, Ri, am at f_adult).
  % An event that cannot be reached at that f (ultimate maturity below the
  % threshold) is returned as NaN; the run itself goes on.
  % get_tx and reprod_rate_foetus freeze when puberty cannot be reached, so
  % (i) they are only called when the event they need is reachable, and
  % (ii) birth and weaning are computed with a puberty target just above
  % weaning, which leaves tau_b and tau_x unchanged (checked numerically:
  % tau_b = 1.92979, l_b = 0.337722 for Sorex at f = 1 with both vectors).
  % (A target just above birth makes get_tx freeze, so it is not used.)
  reach = @(f_val, vH) f_val * max(0, f_val - l_T)^2 > k * vH; % -, DEBtool condition
  pars_tx = [g k l_T v_Hb v_Hx v_Hp];
  pars_x  = [g k l_T v_Hb v_Hx v_Hx*1.001];       % for birth and weaning only

  % life cycle - Piecewise computation for stx mammal model
  if ~reach(f_foetus, v_Hx); info = 0; prdData = []; return; end % f_foetus is not varied in the scenarios
  [~, ~, tau_b, ~, ~, l_b, info] = get_tx(pars_x, f_foetus);
  if ~info; prdData = []; return; end
  tau_x = NaN; tau_bx = NaN; l_x = NaN;
  if reach(f_lactation, v_Hx)
    [~, tau_x, tau_bx, ~, l_x, ~, info_x] = get_tx(pars_x, f_lactation);
    if ~info_x; tau_x = NaN; l_x = NaN; end
  end
  tau_p = NaN; tau_bp = NaN; l_p = NaN;
  if reach(f_juvenile, v_Hp)
    [tau_p, ~, tau_bp, l_p, ~, ~, info_p] = get_tx(pars_tx, f_juvenile);
    if ~info_p; tau_p = NaN; l_p = NaN; end
  end
  % birth
  L_b = L_m * l_b;                           % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_foetus * w);         % g, wet weight at birth at f 
  tT_b = tau_b / k_M / TC;                   % d, time from fertilization to birth
  tT_g = t_0 + tT_b;                         % d, gestation time at f and T

  % weaning
  L_x = L_m * l_x;                           % cm, structural length at weaning
  Ww_x = L_x^3 * (1 + f_lactation * w);      % g, wet weight at weaning
  tT_x = (tau_x - tau_bx)/ k_M / TC;         % d, time since birth at weaning at f and T (NaN if not reached)
  % puberty
  L_p = L_m * l_p;                           % cm, structural length at puberty at f
  Ww_p = L_p^3 *(1 + f_juvenile * w);        % g, wet weight at puberty 
  tT_p = (tau_p - tau_bp)/ k_M / TC;         % d, time since birth at puberty at f and T (NaN if not reached)

  % ultimate
  l_i = f_adult - l_T;                       % -, scaled ultimate length
  L_i = L_m * l_i;                           % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * w);          % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  if reach(f_adult, v_Hp)
    RT_i = TC * reprod_rate_foetus(L_i, f_adult, pars_R);   % #/d, ultimate reproduction rate at T
  else
    RT_i = 0; % maturity maintenance at E_Hp exceeds the flux to maturity/reproduction
  end

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

  % ---------- uni-variate data: ACUTE scenario (n-day window)
  %% --- PIECEWISE FOOD WITH A 7-DAY RESTRICTED STARVATION WINDOW ---
  % Stage onsets = BASELINE events (unrestricted f), so that the window starts
  % at the onset of the stage and the baseline timeline is the same for every
  % reduction factor.
  [~, tau_x0, tau_bx0] = get_tx(pars_x, f_base_lactation);
  [tau_p0, ~, tau_bp0] = get_tx(pars_tx, f_base_juvenile);
  tT_x0 = (tau_x0 - tau_bx0)/ k_M / TC;      % d, baseline time since birth at weaning
  tT_p0 = (tau_p0 - tau_bp0)/ k_M / TC;      % d, baseline time since birth at puberty
  % Full-period scenario (APPENDIX C, par.ct36_full = 1): the reduced f is kept over the
  % whole stage, which then ends at the event reached at the reduced f (Inf if never reached).
  d = 7;                                % d, length of the restriction window (acute scenario)
  if isfield(par, 'ct36_full') && par.ct36_full
    d = Inf;
    if f_lactation ~= f_base_lactation; tT_x0 = tT_x; tT_x0(isnan(tT_x0)) = Inf; end
    if f_juvenile ~= f_base_juvenile; tT_p0 = tT_p; tT_p0(isnan(tT_p0)) = Inf; end
  end
  % 1. Build the unstarved baseline timeline for the whole life cycle
  T = unique([0; tT_x0; tT_p0; min(d, tT_x0); min(tT_x0 + d, tT_p0); tT_p0 + d]); T = T(T >= 0 & isfinite(T)); % d, exact switch times of f
  f_vec = f_base_lactation * ones(size(T));
  f_vec(T >= tT_x0) = f_base_juvenile;
  f_vec(T >= tT_p0) = f_base_adult;
  % 2. Apply the 7-day starvation ONLY if the runner modified the parameter
  if f_lactation ~= f_base_lactation
      f_vec(T >= 0 & T < min(d, tT_x0)) = f_lactation;
  end
  if f_juvenile ~= f_base_juvenile
      f_vec(T >= tT_x0 & T < min(tT_x0 + d, tT_p0)) = f_juvenile;
  end
  if f_adult ~= f_base_adult
      f_vec(T >= tT_p0 & T < tT_p0 + d) = f_adult;
  end
  % Time axis for post-natal development in DEBtool is scaled time since birth
  tau_run = tW(:,1) * k_M * TC;
  % Custom robust ODE integration bypassing get_tx event detection
  el_0 = [1; l_b]; % Reserve density (e) at birth is 1, scaled length is l_b
  el_res = el_path_CT36EFSA(tau_run, T * k_M * TC, f_vec, el_0, g, l_T); % exact n-day windows
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
    % Get current f using step interpolation
    f = interp1(tf_mat(:,1), tf_mat(:,2), tau, 'previous', 'extrap');
    % DEB differential equations for post-natal growth
    de = g * (f - e) / l;
    rho = g * (e/l - 1 - l_T/l) / (e + g);
    dl = l * rho / 3;
    del = [de; dl];
end