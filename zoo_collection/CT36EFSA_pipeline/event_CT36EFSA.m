function [tau_e, tau_b, l_e] = event_CT36EFSA(pars, f, reachable, tel_b)
% event_CT36EFSA: scaled age at a maturity event (metamorphosis/fledging/puberty)
% for the std model at constant food, or NaN if the event cannot be reached.
%
% Input
%   pars:      5-vector [g; k; l_T; v_Hb; v_H_event], as for get_tp
%   f:         scalar, scaled functional response (held constant from birth)
%   reachable: logical, f*(f - l_T)^2 > k*v_H_event (DEBtool condition, get_tp.m);
%              if false, get_tp is not called
%   tel_b:     optional 3-vector [tau_b, e_b, l_b] passed to get_tp (birth state)
%
% Output
%   tau_e: scaled age at the event (NaN if not reached)
%   tau_b: scaled age at birth at the same f (to compute time since birth)
%   l_e:   scaled length at the event (NaN if not reached)
%
% Used by the predict_*_CT36EFSA files so that an unreachable event gives NaN
% for that output only, instead of stopping the whole run.

  tau_e = NaN; tau_b = NaN; l_e = NaN;
  if ~reachable
    return
  end
  if nargin < 4 || isempty(tel_b)
    [tau_e, tau_b, l_e, ~, info] = get_tp(pars, f);
  else
    [tau_e, tau_b, l_e, ~, info] = get_tp(pars, f, tel_b);
  end
  if ~info
    tau_e = NaN; l_e = NaN;
  end
end
