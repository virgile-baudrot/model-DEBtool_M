function run_prediction_Rana_temporaria_CT36EFSA(f_modes_arg)
% run_prediction_Rana_temporaria_CT36EFSA.m
close all; clc;

% 1. Find the absolute path of the current script and parent directory
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(fileparts(current_dir)); % Go up to DEBtool_M
addpath(genpath(parent_dir));

%% 1) Species name
species_name = "Rana_temporaria";

%% 2) Simulation parameters
t_max = 5000; % 14 years lifespan
f_values = 0.5:0.01:1.0; 

if nargin < 1 || isempty(f_modes_arg)
    f_modes = ["f_nest", "f_juvenile", "f_adult"];
else
    f_modes = string(f_modes_arg);
end

%% --- Generate dynamic function names based on species_name
dataFunc    = str2func(sprintf('mydata_%s', species_name));
parsFunc    = str2func(sprintf('pars_init_%s_CT36EFSA', species_name));
predictFunc = str2func(sprintf('predict_%s_CT36EFSA', species_name));

%% 3) Load data and parameters
disp('--> Loading data and parameters...');
[data, auxData, metaData, txtData, weights] = dataFunc();
[par_init, metaPar, txtPar] = parsFunc(metaData);

%% 3b) Base values for each f parameter
f_base.f_nest     = par_init.f_base_nest;
f_base.f_juvenile = par_init.f_base_juvenile;
f_base.f_adult    = par_init.f_base_adult;

%% 4) Prepare the time vector for simulation
times  = linspace(0, t_max, 301)';      
dataSim = data;
dataSim.tL_f = [times, zeros(size(times))]; 
dataSim.tL_m = [times, zeros(size(times))]; 

% Create a smooth theoretical array of Lengths to plot the continuous LdL line
dataSim.LdL = [linspace(2, 9, 100)', zeros(100, 1)];

%% 5) Create a new folder for results
res_dir = fullfile(current_dir, sprintf('results_%s', species_name));
if ~exist(res_dir, 'dir')
    mkdir(res_dir);
end

%% 6) Define color mapping
f_min = min(f_values); f_max = max(f_values); df = 0.01;
f_ref = f_min:df:f_max; nCol = numel(f_ref);
cmap = flipud(turbo(nCol));             
f2idx = @(f) max(1, min(nCol, round((f - f_min)/df) + 1)); 

for mode_idx = 1:numel(f_modes)
    f_mode = f_modes(mode_idx);
    par = par_init; % Reset parameters for this mode
    
    par.f_nest     = f_base.f_nest;
    par.f_juvenile = f_base.f_juvenile;
    par.f_adult    = f_base.f_adult;

    fprintf('\n==================================================\n');
    fprintf('TESTING MODALITY: %s\n', f_mode);
    fprintf('==================================================\n');

    %% 7) Initialize plots
    % --- FEMALE TIME PLOT ---
    fig_f = figure('Name', sprintf('SVL of %s female (%s)', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), ...
                   'Color','w', 'Position', [100, 100, 450, 350]); 
    hold on; box on;
    % No scatter data for time length in Rana temporaria

    % --- MALE TIME PLOT ---
    fig_m = figure('Name', sprintf('SVL of %s male (%s)', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), ...
                   'Color','w', 'Position', [600, 100, 450, 350]); 
    hold on; box on;

    % --- LdL PLOT (Length vs Change in length) ---
    fig_LdL = figure('Name', sprintf('Change in SVL %s (%s)', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), ...
                   'Color','w', 'Position', [1100, 100, 450, 350]); 
    hold on; box on;
    if isfield(data,'LdL') && ~isempty(data.LdL)
        hScatter_LdL = scatter(data.LdL(:,1), data.LdL(:,2), 30, 'k', 'filled', 'DisplayName', 'Data');
    end

    alpha_val = 0.7;      
    lw_default = 2;       
    lw_f1 = 3;            

    vars_to_extract = {'am', 'tj', 'tp', 'Lj', 'Lp', 'Lpm', 'Li', 'Lim', 'Wwb', 'Wwi', 'Wwim', 'Ri'};
    all_results = zeros(numel(f_values), length(vars_to_extract));

    allSim = [];  

    for i = 1:numel(f_values)
        par.(f_mode) = f_base.(f_mode) * f_values(i);
        fprintf('  -> Computing f_value = %.2f...\n', f_values(i));
        
        try
            [prdSim, info] = predictFunc(par, dataSim, auxData);
        catch ME
            fprintf('     [!] ERROR during prediction: %s\n', ME.message);
            info = 0;
        end

        if info ~= 0 && isfield(prdSim,'tL_f') && isfield(prdSim, 'tL_m') && isfield(prdSim, 'LdL')
            idx = f2idx(f_values(i));
            
            % Plot Female SVL
            figure(fig_f);
            if abs(f_values(i) - 1.0) < 1e-6  
                plot(times, prdSim.tL_f, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
            else
                plot(times, prdSim.tL_f, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
            end

            % Plot Male SVL
            figure(fig_m);
            if abs(f_values(i) - 1.0) < 1e-6  
                plot(times, prdSim.tL_m, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
            else
                plot(times, prdSim.tL_m, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
            end

            % Plot Change in SVL (LdL)
            figure(fig_LdL);
            if abs(f_values(i) - 1.0) < 1e-6  
                plot(dataSim.LdL(:,1), prdSim.LdL, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
            else
                plot(dataSim.LdL(:,1), prdSim.LdL, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
            end
            
            for v = 1:length(vars_to_extract)
                var_name = vars_to_extract{v};
                if isfield(prdSim, var_name)
                    all_results(i, v) = prdSim.(var_name);
                else
                    all_results(i, v) = NaN;
                end
            end

            allNames = fieldnames(par);
            names = {}; values = []; units = {}; labels = {};
            for k = 1:numel(allNames)
                nm = allNames{k};
                if strcmp(nm,'free'); continue; end
                val = par.(nm);
                if isnumeric(val) && isscalar(val)
                    names{end+1,1} = nm; values(end+1,1) = val;
                    units{end+1,1} = txtPar.units.(nm);
                    labels{end+1,1} = txtPar.label.(nm);
                end
            end
            paramTable = table(names, values, units, labels, 'VariableNames', {'param','value','unit','label'});
            paramStruct = table2struct(sortrows(paramTable,'param'),'ToScalar',true);

            simStruct = struct('f_value', f_values(i), 'results', prdSim, 'parameters', paramStruct, 'times', times);
            allSim = [allSim; simStruct];          
        else
            fprintf('     [-] Prediction returned info=0. Skipping.\n');
        end
    end

    %% Finalize Female Figure
    figure(fig_f);
    xlabel('Time since birth (d)', 'Color', 'k'); ylabel('SVL (cm)', 'Color', 'k');
    title(sprintf('%s female SVL - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color', 'k');
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',11); grid on;
    colormap(cmap); cb = colorbar; cb.Label.String = strrep(f_mode,'_','\_'); cb.Color = 'k'; caxis([f_min f_max]);
    
    figFilenamePNG_f = fullfile(res_dir, sprintf('%s_%s_female.png', species_name, f_mode));
    figFilenameJPG_f = fullfile(res_dir, sprintf('%s_%s_female.jpg', species_name, f_mode));
    exportgraphics(fig_f, figFilenamePNG_f, 'BackgroundColor', 'white', 'Resolution', 150);
    exportgraphics(fig_f, figFilenameJPG_f, 'BackgroundColor', 'white', 'Resolution', 150);

    %% Finalize Male Figure
    figure(fig_m);
    xlabel('Time since birth (d)', 'Color', 'k'); ylabel('SVL (cm)', 'Color', 'k');
    title(sprintf('%s male SVL - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color', 'k');
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',11); grid on;
    colormap(cmap); cb = colorbar; cb.Label.String = strrep(f_mode,'_','\_'); cb.Color = 'k'; caxis([f_min f_max]);
    
    figFilenamePNG_m = fullfile(res_dir, sprintf('%s_%s_male.png', species_name, f_mode));
    figFilenameJPG_m = fullfile(res_dir, sprintf('%s_%s_male.jpg', species_name, f_mode));
    exportgraphics(fig_m, figFilenamePNG_m, 'BackgroundColor', 'white', 'Resolution', 150);
    exportgraphics(fig_m, figFilenameJPG_m, 'BackgroundColor', 'white', 'Resolution', 150);

    %% Finalize LdL Figure
    figure(fig_LdL);
    if exist('hScatter_LdL', 'var'), uistack(hScatter_LdL,'top'); end
    xlabel('SVL (cm)', 'Color', 'k'); ylabel('Change in SVL (cm/d)', 'Color', 'k');
    title(sprintf('%s LdL - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color', 'k');
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',11); grid on;
    colormap(cmap); cb = colorbar; cb.Label.String = strrep(f_mode,'_','\_'); cb.Color = 'k'; caxis([f_min f_max]);
    
    % --- Crop the X-axis (keep minimum as is, set maximum to 6.5) ---
    xlim([-inf, 6.5]);
    
    figFilenamePNG_LdL = fullfile(res_dir, sprintf('%s_%s_LdL.png', species_name, f_mode));
    figFilenameJPG_LdL = fullfile(res_dir, sprintf('%s_%s_LdL.jpg', species_name, f_mode));
    exportgraphics(fig_LdL, figFilenamePNG_LdL, 'BackgroundColor', 'white', 'Resolution', 150);
    exportgraphics(fig_LdL, figFilenameJPG_LdL, 'BackgroundColor', 'white', 'Resolution', 150);

    %% 8) Save all simulations
    jsonStr = jsonencode(allSim);
    jsonFilename = fullfile(res_dir, sprintf('simulations_CT36EFSA_%s_%s.json', species_name, f_mode));
    fid = fopen(jsonFilename,'w'); fwrite(fid, jsonStr, 'char'); fclose(fid);

    %% 9) Save extracted outputs
    results_table = array2table(all_results, 'VariableNames', vars_to_extract);
    results_table.f_value = f_values';  
    results_table = [results_table(:, end), results_table(:, 1:end-1)];
    writetable(results_table, fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_mode)));

end

end