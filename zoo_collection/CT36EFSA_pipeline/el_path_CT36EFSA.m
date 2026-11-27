function el = el_path_CT36EFSA(tau_out, tau_sw, f_sw, el0, g, l_T)
% el_path_CT36EFSA: scaled reserve e and length l along a piecewise-constant food path (CT36 scenarios)
%
% tau_out : scaled times (ascending, first = time of el0) at which e and l are returned
% tau_sw  : scaled times at which f changes (ascending; the first one <= tau_out(1))
% f_sw    : f from each switch time onward
% el0     : [e; l] at tau_out(1)
%
% The ODE is integrated segment by segment between the exact switch times, so a
% restriction window keeps its exact length whatever the output grid (sampling f on
% the output grid shortened, lengthened or skipped the n-day windows).
% Same equations as the former dget_el_local (std model, structure may shrink).

  tau_out = tau_out(:); tau_sw = tau_sw(:); f_sw = f_sw(:);
  el = NaN(numel(tau_out), 2);
  opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
  y = el0(:);
  edges = [tau_sw; max(tau_out(end), tau_sw(end)) + 1];
  for s = 1:numel(tau_sw)
    a = max(edges(s), tau_out(1)); b = edges(s+1);
    if b <= a; continue; end
    if s == numel(tau_sw); idx = find(tau_out >= a & tau_out <= b); else; idx = find(tau_out >= a & tau_out < b); end
    pts = unique([a; tau_out(idx); b]);
    f = f_sw(s);
    rhs = @(t, x) del_el(x, f, g, l_T);
    if numel(pts) == 2
      [~, Y] = ode45(rhs, pts, y, opts); Y = Y([1 end], :);
    else
      [~, Y] = ode45(rhs, pts, y, opts);
    end
    [~, loc] = ismember(tau_out(idx), pts);
    el(idx, :) = Y(loc, :);
    y = Y(end, :)';
  end
end

function d = del_el(x, f, g, l_T)
  e = max(1e-6, x(1)); l = max(1e-6, x(2));
  de = g * (f - e) / l;
  rho = g * (e/l - 1 - l_T/l) / (e + g);
  d = [de; l * rho / 3];
end
