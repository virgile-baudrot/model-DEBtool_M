function [x, fval, info, output] = fzero(fun, x0, options, varargin)
  if nargin < 3 || isempty(options), options = optimset(); end
  if ~isempty(varargin), f = @(x) fun(x, varargin{:}); else, f = fun; end
  [x, fval, info, output] = fzero_orig(f, x0, options);
end
