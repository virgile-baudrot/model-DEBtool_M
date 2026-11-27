% run_octave_CT36EFSA.m -- CT36 EFSA food-reduction simulations (GNU Octave 8.4).
% For each stage, f_stage = f_base_stage * factor, factor = 0.5:0.01:1.0, other stages
% at their baseline. Two scenarios:
%   acute : the reduced f lasts an n-day window from the onset of the stage
%           (report section 3.3, top panels of Figs 3.2-3.11)
%   full  : the reduced f lasts the whole stage (report APPENDIX C, Figs C.1-C.12);
%           sets par.ct36_full = 1 in the predict_*_CT36EFSA files
% The chronic outputs (outputs_*.csv, bottom panels of Figs 3.3-3.11) do not depend on
% the scenario and are written by the acute run only.
%
% Writes, in <outroot>/<species>/:
%   outputs_<species>_<mode>.csv  (acute only) same columns as the MATLAB CSV (f_value, then variables;
%                                 empty cell = not computed / event not reached)
%   traj_<species>_<mode>.csv     (acute; traj_<species>_<mode>_full.csv for full) long format: f_value,panel,x,y  (growth curves)
%   data_<species>.csv            observed data drawn on the growth panels: panel,x,y,source
% Figures are drawn afterwards by plot_CT36EFSA.R.
%
% Usage (shell):  octave run_octave_CT36EFSA.m <species> [outroot] [acute|full]
% (run_all_CT36EFSA.sh runs everything)

crash_dumps_octave_core(false);
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));                       % DEBtool_M
addpath(genpath(root));
rmpath(fullfile(here, 'octave_shims'));
if exist('OCTAVE_VERSION', 'builtin'); addpath(fullfile(here, 'octave_shims'), '-begin'); end % Octave only
warning('off', 'all');

args = argv();
sp = args{1};
if numel(args) >= 2 && ~isempty(args{2}); outroot = args{2}; else; outroot = fullfile(root, 'zoo_collection', 'CT36EFSA_results'); end
if numel(args) >= 3; scenario = args{3}; else; scenario = 'acute'; end
if ~any(strcmp(scenario, {'acute', 'full'})); error('scenario must be acute or full'); end
is_full = strcmp(scenario, 'full'); sfx = ''; if is_full; sfx = '_full'; end
outdir = fullfile(outroot, sp);
if ~exist(outdir, 'dir'); mkdir(outdir); end

cd(fullfile(root, 'zoo_collection', sp));
[data, auxData, metaData] = feval(['mydata_' sp]);
par_init = feval(['pars_init_' sp '_CT36EFSA'], metaData);
par_init.ct36_full = double(is_full);
dataSim = data;

% --- species settings (copied from the MATLAB runners)
switch sp
  case {'Perdix_perdix', 'Motacilla_flava', 'Alauda_arvensis', 'Turdus_philomelos'}
    modes = {'f_nest', 'f_juvenile', 'f_adult'};
    times = linspace(0, 150, 301)'; dataSim.tW = [times, 0*times];
    if strcmp(sp, 'Perdix_perdix'); vars = {'ab','tx','tp','tR','am','Ww0','Wwb','Wwi','Wwim','Ri'};
    else; vars = {'ab','tx','tp','tR','am','Wwb','Wwi','Ri'}; end
  case 'Sorex_araneus'
    modes = {'f_lactation', 'f_juvenile', 'f_adult'};
    times = linspace(0, 500, 301)'; dataSim.tW = [times, 0*times];
    vars = {'tg','tx','tp','am','Li','Wwb','Wwx','Wwp','Wwi','Ri'};
  case 'Crocidura_russula'
    modes = {'f_lactation', 'f_juvenile', 'f_adult'};
    times = linspace(0, 500, 301)'; dataSim.tW = [times, 0*times];
    vars = {'tg','tx','tp','am','Wwb','Wwx','Wwp','Wwi','Ri'};
  case 'Ichthyosaura_alpestris'
    modes = {'f_nest', 'f_juvenile', 'f_adult'};
    times = linspace(0, 8000, 301)'; dataSim.tL_f = [times, 0*times]; dataSim.tL_m = [times, 0*times];
    vars = {'am','tj','tp','Lj','Lp','Lpm','Li','Lim','Wwb','Wwi','Wwim','Ri'};
  case 'Rana_temporaria'
    modes = {'f_nest', 'f_juvenile', 'f_adult'};
    times = linspace(0, 5000, 301)'; dataSim.tL_f = [times, 0*times]; dataSim.tL_m = [times, 0*times];
    dataSim.LdL = [linspace(2, 9, 100)', zeros(100, 1)];
    vars = {'am','tj','tp','Lj','Lp','Lpm','Li','Lim','Wwb','Wwi','Wwim','Ri'};
  case 'Podarcis_muralis'
    modes = {'f_embryo', 'f_juvenile', 'f_adult'};
    dataSim.t_run = linspace(0, 6000, 300)'; dataSim.L_run = linspace(2, 7.5, 300)';
    vars = {'tp','am','Lb','Lp','Lpm','Li','Lim','Ww0','Wwb','Wwi','Wwim','Ri'};
  otherwise
    error('unknown species %s', sp);
end

% --- observed data for the growth panels
if ~is_full % observed data are the same for both scenarios
fid = fopen(fullfile(outdir, sprintf('data_%s.csv', sp)), 'w');
fprintf(fid, 'panel,x,y,source\n');
put = @(fid, panel, x, y, src) arrayfun(@(i) fprintf(fid, '%s,%.10g,%.10g,%s\n', panel, x(i), y(i), src), 1:numel(x));
if isfield(data, 'tW');   put(fid, 'tW',   data.tW(:,1),   data.tW(:,2),   'data'); end
if isfield(data, 'tL_f'); put(fid, 'tL_f', data.tL_f(:,1), data.tL_f(:,2), 'data'); end
if isfield(data, 'tL_m'); put(fid, 'tL_m', data.tL_m(:,1), data.tL_m(:,2), 'data'); end
if isfield(data, 'LdL');  put(fid, 'LdL',  data.LdL(:,1),  data.LdL(:,2),  'data'); end
if isfield(data, 'tL_K'); put(fid, 'tL_f', data.tL_K(:,1), data.tL_K(:,2), 'tL_K');
  if size(data.tL_K, 2) >= 3; put(fid, 'tL_m', data.tL_K(:,1), data.tL_K(:,3), 'tL_K'); end; end
if isfield(data, 'tL_D'); put(fid, 'tL_f', data.tL_D(:,1), data.tL_D(:,2), 'tL_D');
  if size(data.tL_D, 2) >= 3; put(fid, 'tL_m', data.tL_D(:,1), data.tL_D(:,3), 'tL_D'); end; end
if isfield(data, 'LW_f'); put(fid, 'LW_f', data.LW_f(:,1), data.LW_f(:,2), 'data'); end
if isfield(data, 'LW_m'); put(fid, 'LW_m', data.LW_m(:,1), data.LW_m(:,2), 'data'); end
if isfield(data, 'LN');   put(fid, 'LN',   data.LN(:,1),   data.LN(:,2),   'data'); end
fclose(fid);
end

f_values = 0.5:0.01:1.0;
for m = 1:numel(modes)
  md = modes{m}; par = par_init;
  % all stages at baseline, then the active one is scaled
  stages = fieldnames(par_init); stages = stages(strncmp(stages, 'f_base_', 7));
  for s = 1:numel(stages); par.(strrep(stages{s}, 'f_base_', 'f_')) = par_init.(stages{s}); end
  base = par_init.(['f_base_' md(3:end)]);
  res = NaN(numel(f_values), numel(vars));
  ftr = fopen(fullfile(outdir, sprintf('traj_%s_%s%s.csv', sp, md, sfx)), 'w');
  fprintf(ftr, 'f_value,panel,x,y\n');
  t0 = tic; nok = 0;
  for i = 1:numel(f_values)
    par.(md) = base * f_values(i);
    try
      [prd, info] = feval(['predict_' sp '_CT36EFSA'], par, dataSim, auxData);
    catch err
      printf('  ERROR %s %.2f: %s\n', md, f_values(i), err.message); info = 0;
    end
    if ~info; continue; end
    nok = nok + 1;
    for v = 1:numel(vars); if isfield(prd, vars{v}) && ~isempty(prd.(vars{v})); res(i,v) = prd.(vars{v}); end; end
    % growth curves
    P = {};
    if isfield(prd, 'tW');   P(end+1,:) = {'tW',   times, prd.tW};   end
    if isfield(prd, 'tL_f'); P(end+1,:) = {'tL_f', times, prd.tL_f}; end
    if isfield(prd, 'tL_m'); P(end+1,:) = {'tL_m', times, prd.tL_m}; end
    if isfield(prd, 'LdL');  P(end+1,:) = {'LdL',  dataSim.LdL(:,1), prd.LdL}; end
    if isfield(prd, 'plot_t')
      P(end+1,:) = {'tL_f', prd.plot_t,     prd.plot_SVL_f};
      P(end+1,:) = {'tL_m', prd.plot_t,     prd.plot_SVL_m};
      P(end+1,:) = {'LW_f', prd.plot_SVL_f, prd.plot_Ww_f};
      P(end+1,:) = {'LW_m', prd.plot_SVL_m, prd.plot_Ww_m};
      P(end+1,:) = {'LN',   prd.plot_L_run, prd.plot_LN};
    end
    for p = 1:size(P, 1)
      x = P{p,2}(:); y = P{p,3}(:);
      for k = 1:numel(x); fprintf(ftr, '%.2f,%s,%.6g,%.6g\n', f_values(i), P{p,1}, x(k), y(k)); end
    end
  end
  fclose(ftr);
  if ~is_full
  fid = fopen(fullfile(outdir, sprintf('outputs_%s_%s.csv', sp, md)), 'w');
  fprintf(fid, 'f_value,%s\n', strjoin(vars, ','));
  for i = 1:numel(f_values)
    fprintf(fid, '%.2f', f_values(i));
    for v = 1:numel(vars)
      if isnan(res(i,v)); fprintf(fid, ','); else; fprintf(fid, ',%.10g', res(i,v)); end
    end
    fprintf(fid, '\n');
  end
  fclose(fid);
  end
  printf('%s %s %s (base %.4g): %d/%d runs ok, %.1fs\n', sp, scenario, md, base, nok, numel(f_values), toc(t0)); fflush(stdout);
end
