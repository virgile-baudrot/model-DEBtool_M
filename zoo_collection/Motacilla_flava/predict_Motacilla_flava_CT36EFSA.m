function [prdData, info] = predict_Motacilla_flava_CT36EFSA(par, data, auxData)

  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par);
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC = tempcorr(temp.am, T_ref, T_A); kT_M = TC * k_M;

  % ---------- zero-variate data: CHRONIC outputs
  % Each stage f is held constant from birth (tx at f_nest, tp at f_juvenile)
  % or over the whole life (Li, Wwi, Ri, am at f_adult).
  % An event that cannot be reached at that f (ultimate maturity below the
  % threshold, see get_tp.m) is returned as NaN; the run itself goes on.
  reach = @(f_val, vH) f_val * max(0, f_val - l_T)^2 > k * vH; % -, DEBtool condition
  pars_tx = [g; k; l_T; v_Hb; v_Hx];  % compose parameter vector
  pars_tp = [g; k; l_T; v_Hb; v_Hp];  % compose parameter vector
  [tau_b, l_b, info] = get_tb([g; k; v_Hb], f_embryo); % -, scaled age & length at birth
  if ~info; prdData = []; return; end
  tau_x = NaN; tau_bx = NaN; l_x = NaN;
  if reach(f_nest, v_Hx)
    [tau_x, tau_bx, l_x, ~, info_x] = get_tp(pars_tx, f_nest);
    if ~info_x; tau_x = NaN; l_x = NaN; end
  end
  tau_p = NaN; tau_bp = NaN; l_p = NaN;
  if reach(f_juvenile, v_Hp)
    [tau_p, tau_bp, l_p, ~, info_p] = get_tp(pars_tp, f_juvenile);
    if ~info_p; tau_p = NaN; l_p = NaN; end
  end

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * (1 + f_embryo * ome);       % g, wet weight at birth at f
  tT_b = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T
  aT_b = tT_b;                % d, time at birth at f and T
  t_g = tT_b;                 % d, gestation time (time from fertilization to birth)
  % fledging
  tT_x = (tau_x - tau_bx)/ kT_M;        % d, time since birth at fledging (NaN if not reached)
  L_x = L_m * l_x;
  Ww_x = L_x^3 * (1 + f_nest * ome);

  % puberty
  tT_p = (tau_p - tau_bp)/ kT_M;        % d, time since birth at puberty (NaN if not reached)
  L_p = L_m * l_p;
  Ww_p = L_p^3 * (1 + f_juvenile * ome);

  % ultimate
  l_i = f_adult - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Ww_i = L_i^3 * (1 + f_adult * ome);       % g, ultimate wet weight
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  if reach(f_adult, v_Hp)
    RT_i = TC * reprod_rate(L_i, f_adult, pars_R);                % #/d, ultimate reproduction rate at T
  else
    RT_i = 0; % maturity maintenance at E_Hp exceeds the flux to maturity/reproduction
  end

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

  % ---------- uni-variate data: ACUTE scenario (n-day window)

  %% --- PIECEWISE FOOD WITH A 7-DAY RESTRICTED STARVATION WINDOW ---
  % Stage onsets = BASELINE events (unrestricted f), so that the window starts
  % at the onset of the stage and the baseline timeline is the same for every
  % reduction factor.
  [tau_x0, tau_bx0] = get_tp(pars_tx, f_base_nest);
  [tau_p0, tau_bp0] = get_tp(pars_tp, f_base_juvenile);
  tT_x0 = (tau_x0 - tau_bx0)/ kT_M;     % d, baseline time since birth at fledging
  tT_p0 = (tau_p0 - tau_bp0)/ kT_M;     % d, baseline time since birth at puberty
  % Full-period scenario (APPENDIX C, par.ct36_full = 1): the reduced f is kept over the
  % whole stage, which then ends at the event reached at the reduced f (Inf if never reached).
  d = 7;                                % d, length of the restriction window (acute scenario)
  if isfield(par, 'ct36_full') && par.ct36_full
    d = Inf;
    if f_nest ~= f_base_nest; tT_x0 = tT_x; tT_x0(isnan(tT_x0)) = Inf; end
    if f_juvenile ~= f_base_juvenile; tT_p0 = tT_p; tT_p0(isnan(tT_p0)) = Inf; end
  end

  % 1. Build the unstarved baseline timeline for the whole life cycle
  f = f_base_nest * ones(size(tW(:,1)));
  f(tW(:,1) >= tT_x0) = f_base_juvenile;
  f(tW(:,1) >= tT_p0) = f_base_adult;

  % 2. Apply the 7-day starvation ONLY if the runner modified the parameter
  % Nest stage starvation (starts at birth t=0, lasts 7 days or until fledging)
  if f_nest ~= f_base_nest
      f(tW(:,1) >= 0 & tW(:,1) < min(d, tT_x0)) = f_nest;
  end
  % Juvenile stage starvation (starts at fledging, lasts 7 days)
  if f_juvenile ~= f_base_juvenile
      f(tW(:,1) >= tT_x0 & tW(:,1) < min(tT_x0 + d, tT_p0)) = f_juvenile;
  end
  % Adult stage starvation (starts at puberty, lasts 7 days)
  if f_adult ~= f_base_adult
      f(tW(:,1) >= tT_p0 & tW(:,1) < tT_p0 + d) = f_adult;
  end

  % Build time-food table expected by get_tp
  tf = [tW(:,1) * kT_M, f];

  % Call get_tp once (vectorized in time, piecewise in f)
  tvel = get_tp(pars_tp, tf, [tau_b 1 l_b], tW(:,1) * kT_M);
  % Compute wet weight
  EWw = (L_m .* tvel(:,4)).^3 .* (1 + tvel(:,3) .* ome);

  % pack to output
  prdData.tW = EWw;

end
