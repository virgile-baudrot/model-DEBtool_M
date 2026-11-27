function [q, n] = quad(fun, a, b, tol, trace, varargin)
  if nargin < 4 || isempty(tol), tol = 1e-6; end
  g = @(x) arrayfun(@(xx) fun(xx, varargin{:}), x);
  q = integral(g, a, b, 'AbsTol', tol); n = 0;
end
