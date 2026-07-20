function run_prediction_Podarcis_muralis_CT36EFSA(f_modes_arg)
% run_prediction_Podarcis_muralis_CT36EFSA.m
close all; clc;

% 1. Find the absolute path of the current script and parent directory
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(fileparts(current_dir)); 
addpath(genpath(parent_dir));

%% 1) Species name
species_name = "Podarcis_muralis";

%% 2) Simulation parameters
t_max = 6000; % Time array (days)
L_max = 7.5;  % Length array (cm)
f_values = 0.5:0.01:1.0; 

if nargin < 1 || isempty(f_modes_arg)
    f_modes = ["f_embryo", "f_juvenile", "f_adult"];
else
    f_modes = string(f_modes_arg);
end

%% --- Generate dynamic function names
dataFunc    = str2func(sprintf('mydata_%s', species_name));
parsFunc    = str2func(sprintf('pars_init_%s_CT36EFSA', species_name));
predictFunc = str2func(sprintf('predict_%s_CT36EFSA', species_name));

%% 3) Load data and parameters
disp('--> Loading data and parameters...');
[data, auxData, metaData, txtData, weights] = dataFunc();
[par_init, metaPar, txtPar] = parsFunc(metaData);

f_base.f_embryo   = par_init.f_base_embryo;
f_base.f_juvenile = par_init.f_base_juvenile;
f_base.f_adult    = par_init.f_base_adult;

%% 4) Prepare continuous vectors for simulation plots
dataSim = data;
% Provide the boundaries for continuous plot integration inside predict
dataSim.t_run = linspace(0, t_max, 300)';      
dataSim.L_run = linspace(2, L_max, 300)';

%% 5) Create a new folder for results
res_dir = fullfile(current_dir, sprintf('results_%s', species_name));
if ~exist(res_dir, 'dir'), mkdir(res_dir); end

%% 6) Define color mapping
f_min = min(f_values); f_max = max(f_values); df = 0.01;
f_ref = f_min:df:f_max; nCol = numel(f_ref);
cmap = flipud(turbo(nCol));             
f2idx = @(f) max(1, min(nCol, round((f - f_min)/df) + 1)); 

for mode_idx = 1:numel(f_modes)
    f_mode = f_modes(mode_idx);
    par = par_init; 
    
    par.f_embryo   = f_base.f_embryo;
    par.f_juvenile = f_base.f_juvenile;
    par.f_adult    = f_base.f_adult;

    fprintf('\n==================================================\n');
    fprintf('TESTING MODALITY: %s\n', f_mode);
    fprintf('==================================================\n');

    %% 7) Initialize plots (Height 220px to fit well in grid)
    pos_w = 400; pos_h = 220;
    
    % SVL vs Time (Female)
    fig_tL_f = figure('Name', 'SVL vs Time (Female)', 'Color','w', 'Position', [100, 700, pos_w, pos_h]); hold on; box on;
    if isfield(data,'tL_K') && ~isempty(data.tL_K), hScat_tLK_f = scatter(data.tL_K(:,1), data.tL_K(:,2), 20, 'k', 'filled'); end
    if isfield(data,'tL_D') && ~isempty(data.tL_D), hScat_tLD_f = scatter(data.tL_D(:,1), data.tL_D(:,2), 20, 'b', 'filled'); end
    
    % SVL vs Time (Male)
    fig_tL_m = figure('Name', 'SVL vs Time (Male)', 'Color','w', 'Position', [550, 700, pos_w, pos_h]); hold on; box on;
    if isfield(data,'tL_K') && size(data.tL_K, 2) >= 3, hScat_tLK_m = scatter(data.tL_K(:,1), data.tL_K(:,3), 20, 'k', 'filled'); end
    if isfield(data,'tL_D') && size(data.tL_D, 2) >= 3, hScat_tLD_m = scatter(data.tL_D(:,1), data.tL_D(:,3), 20, 'b', 'filled'); end

    % Weight vs SVL (Female)
    fig_LW_f = figure('Name', 'Weight vs SVL (Female)', 'Color','w', 'Position', [100, 400, pos_w, pos_h]); hold on; box on;
    if isfield(data,'LW_f') && ~isempty(data.LW_f), hScat_LW_f = scatter(data.LW_f(:,1), data.LW_f(:,2), 20, 'k', 'filled'); end

    % Weight vs SVL (Male)
    fig_LW_m = figure('Name', 'Weight vs SVL (Male)', 'Color','w', 'Position', [550, 400, pos_w, pos_h]); hold on; box on;
    if isfield(data,'LW_m') && ~isempty(data.LW_m), hScat_LW_m = scatter(data.LW_m(:,1), data.LW_m(:,2), 20, 'k', 'filled'); end

    % Fecundity vs SVL (Female)
    fig_LN = figure('Name', 'Fecundity vs SVL (Female)', 'Color','w', 'Position', [1000, 400, pos_w, pos_h]); hold on; box on;
    if isfield(data,'LN') && ~isempty(data.LN), hScat_LN = scatter(data.LN(:,1), data.LN(:,2), 20, 'k', 'filled'); end

    alpha_val = 0.7; lw_def = 2; lw_f1 = 3;            
    vars_to_ext = {'tp', 'am', 'Lb', 'Lp', 'Lpm', 'Li', 'Lim', 'Ww0', 'Wwb', 'Wwi', 'Wwim', 'Ri'};
    all_results = zeros(numel(f_values), length(vars_to_ext));
    allSim = [];  

    for i = 1:numel(f_values)
        par.(f_mode) = f_base.(f_mode) * f_values(i);
        fprintf('  -> Computing f_value = %.2f...\n', f_values(i));
        
        try
            [prdSim, info] = predictFunc(par, dataSim, auxData);
        catch ME
            fprintf('     [!] ERROR: %s\n', ME.message); info = 0;
        end

        if info ~= 0 && isfield(prdSim,'plot_t')
            idx = f2idx(f_values(i));
            lw = lw_def; col = [cmap(idx,:) alpha_val];
            if abs(f_values(i) - 1.0) < 1e-6, lw = lw_f1; col = cmap(idx,:); end
            
            % Plot continuous lines
            figure(fig_tL_f); plot(prdSim.plot_t, prdSim.plot_SVL_f, 'LineWidth', lw, 'Color', col, 'HandleVisibility','off');
            figure(fig_tL_m); plot(prdSim.plot_t, prdSim.plot_SVL_m, 'LineWidth', lw, 'Color', col, 'HandleVisibility','off');
            figure(fig_LW_f); plot(prdSim.plot_SVL_f, prdSim.plot_Ww_f, 'LineWidth', lw, 'Color', col, 'HandleVisibility','off');
            figure(fig_LW_m); plot(prdSim.plot_SVL_m, prdSim.plot_Ww_m, 'LineWidth', lw, 'Color', col, 'HandleVisibility','off');
            figure(fig_LN);   plot(prdSim.plot_L_run, prdSim.plot_LN, 'LineWidth', lw, 'Color', col, 'HandleVisibility','off');
            
            % Extract variables for JSON/CSV
            for v = 1:length(vars_to_ext)
                if isfield(prdSim, vars_to_ext{v}), all_results(i, v) = prdSim.(vars_to_ext{v}); else, all_results(i, v) = NaN; end
            end

            allNames = fieldnames(par); names = {}; values = []; units = {}; labels = {};
            for k = 1:numel(allNames)
                nm = allNames{k};
                if strcmp(nm,'free'); continue; end
                if isnumeric(par.(nm)) && isscalar(par.(nm))
                    names{end+1,1} = nm; values(end+1,1) = par.(nm);
                    units{end+1,1} = txtPar.units.(nm); labels{end+1,1} = txtPar.label.(nm);
                end
            end
            paramStruct = table2struct(sortrows(table(names, values, units, labels, 'VariableNames', {'param','value','unit','label'}),'param'),'ToScalar',true);
            allSim = [allSim; struct('f_value', f_values(i), 'results', prdSim, 'parameters', paramStruct)];          
        else
            fprintf('     [-] Prediction returned info=0. Skipping.\n');
        end
    end

    %% Finalize and Save all 5 plots
    plots_info = {
        fig_tL_f, 'Time since birth (d)', 'SVL (cm)', 'SVL vs Time (Female)', 'tL_female';
        fig_tL_m, 'Time since birth (d)', 'SVL (cm)', 'SVL vs Time (Male)', 'tL_male';
        fig_LW_f, 'SVL (cm)', 'Wet Weight (g)', 'Weight vs SVL (Female)', 'LW_female';
        fig_LW_m, 'SVL (cm)', 'Wet Weight (g)', 'Weight vs SVL (Male)', 'LW_male';
        fig_LN,   'SVL (cm)', 'Yearly Fecundity (#)', 'Fecundity vs SVL (Female)', 'LN';
    };

    for p = 1:size(plots_info, 1)
        figure(plots_info{p, 1});
        xlabel(plots_info{p, 2}, 'Color', 'k'); ylabel(plots_info{p, 3}, 'Color', 'k');
        title(sprintf('%s - %s', plots_info{p, 4}, strrep(f_mode,'_',' ')), 'Color', 'k');
        set(gca,'Color','w','XColor','k','YColor','k','FontSize',10); grid on;
        colormap(cmap); cb = colorbar; cb.Label.String = strrep(f_mode,'_','\_'); cb.Color = 'k'; caxis([f_min f_max]);
        
        % EXPORT GRAPHICS (Robust to dark mode and indexing)
        exportgraphics(plots_info{p, 1}, fullfile(res_dir, sprintf('%s_%s_%s.png', species_name, f_mode, plots_info{p, 5})), 'BackgroundColor', 'white', 'Resolution', 150);
        exportgraphics(plots_info{p, 1}, fullfile(res_dir, sprintf('%s_%s_%s.jpg', species_name, f_mode, plots_info{p, 5})), 'BackgroundColor', 'white', 'Resolution', 150);
    end

    %% Save Data
    fid = fopen(fullfile(res_dir, sprintf('simulations_CT36EFSA_%s_%s.json', species_name, f_mode)),'w');
    fwrite(fid, jsonencode(allSim), 'char'); fclose(fid);

    results_table = array2table(all_results, 'VariableNames', vars_to_ext);
    results_table.f_value = f_values';  
    results_table = [results_table(:, end), results_table(:, 1:end-1)];
    writetable(results_table, fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_mode)));
end
end