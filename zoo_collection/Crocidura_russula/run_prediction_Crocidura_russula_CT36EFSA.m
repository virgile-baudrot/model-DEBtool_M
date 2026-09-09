function run_prediction_Crocidura_russula_CT36EFSA(f_modes_arg)
% run_prediction_Crocidura_russula_CT36EFSA.m
close all; clc;

% 1. Find the absolute path of the current script and parent directory
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(fileparts(current_dir)); % Go up to DEBtool_M
addpath(genpath(parent_dir));

%% 1) Species name
species_name = "Crocidura_russula";

%% 2) Simulation parameters
t_max = 500; % Shrews are short lived, 500 days is a good upper bound
% f_values = reducing factor in [0.5, 1.0] applied to base values below
f_values = 0.5:0.01:1.0;

% --- Which f is varied in the simulations (Mammals logic)
if nargin < 1 || isempty(f_modes_arg)
    f_modes = ["f_lactation", "f_juvenile", "f_adult"];
else
    f_modes = string(f_modes_arg);
end

%% --- Generate dynamic function names based on species_name
dataFunc    = str2func(sprintf('mydata_%s', species_name));
parsFunc    = str2func(sprintf('pars_init_%s_CT36EFSA', species_name));
predictFunc = str2func(sprintf('predict_%s_CT36EFSA', species_name));

%% 3) Load data and parameters
[data, auxData, metaData, txtData, weights] = dataFunc();
[par_init, metaPar, txtPar] = parsFunc(metaData);

%% 3b) Base values for each f parameter (reducing factor applies to these)
f_base.f_lactation = par_init.f_base_lactation;
f_base.f_juvenile  = par_init.f_base_juvenile;
f_base.f_adult     = par_init.f_base_adult;

%% 4) Prepare the time vector for simulation
times  = linspace(0, t_max, 301)';      % 301 equally spaced points
dataSim = data;
dataSim.tW = [times, zeros(size(times))]; % placeholder for simulated weights

%% 5) Create a new folder for results
res_dir = fullfile(current_dir, sprintf('results_%s', species_name));
if ~exist(res_dir, 'dir')
    mkdir(res_dir);
end

%% 6) Define color mapping
f_min = min(f_values); f_max = max(f_values); df = 0.01;
f_ref = f_min:df:f_max; nCol = numel(f_ref);
cmap = flipud(turbo(nCol));             % flip turbo so blue = high f
f2idx = @(f) max(1, min(nCol, round((f - f_min)/df) + 1)); 


for mode_idx = 1:numel(f_modes)
    f_mode = f_modes(mode_idx);
    par = par_init; % Reset parameters for this mode
    
    % Lock all f parameters at their base values (non-active ones stay fixed)
    par.f_lactation = f_base.f_lactation;
    par.f_juvenile  = f_base.f_juvenile;
    par.f_adult     = f_base.f_adult;

%% 7) Plot predictions + experimental data
    figure( ...
        'Name', sprintf('Growth of %s (%s)', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), ...
        'Color','w', ...
        'Position', [100, 100, 450, 350]); 
    hold on; box on;

    % --- Scatter experimental data on top
    if isfield(data,'tW') && ~isempty(data.tW)
        hScatter = scatter(data.tW(:,1), data.tW(:,2), 30, 'k', 'filled', 'DisplayName', 'Data');
    end

    % --- Loop over f values and plot predictions & extract variables
    alpha_val = 0.7;      
    lw_default = 2;       
    lw_f1 = 3;            

    % Initialize array for the table based on mammal mapping
    vars_to_extract = {'tg', 'tx', 'tp', 'am', 'Wwb', 'Wwx', 'Wwp', 'Wwi', 'Ri'};
    all_results = zeros(numel(f_values), length(vars_to_extract));

    allSim = [];  % initialize struct array for JSON

    for i = 1:numel(f_values)
        % Apply reducing factor to the base value of the active f parameter
        par.(f_mode) = f_base.(f_mode) * f_values(i);
        [~, prdSim, info] = evalc('predictFunc(par, dataSim, auxData)');

        if info ~= 0 && isfield(prdSim,'tW')
            idx = f2idx(f_values(i));
            if abs(f_values(i) - 1.0) < 1e-6  
                plot(times, prdSim.tW, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
            else
                plot(times, prdSim.tW, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
            end
            
            % Extract variables for CSV
            for v = 1:length(vars_to_extract)
                var_name = vars_to_extract{v};
                if isfield(prdSim, var_name)
                    all_results(i, v) = prdSim.(var_name);
                else
                    all_results(i, v) = NaN;
                end
            end

            % Construct parameter table for JSON
            allNames = fieldnames(par);
            names = {}; values = []; units = {}; labels = {};
            for k = 1:numel(allNames)
                nm = allNames{k};
                if strcmp(nm,'free'); continue; end
                val = par.(nm);
                if isnumeric(val) && isscalar(val)
                    names{end+1,1} = nm;
                    values(end+1,1) = val;
                    if isfield(txtPar,'units') && isfield(txtPar.units,nm)
                        units{end+1,1} = txtPar.units.(nm);
                    else
                        units{end+1,1} = '';
                    end
                    if isfield(txtPar,'label') && isfield(txtPar.label,nm)
                        labels{end+1,1} = txtPar.label.(nm);
                    else
                        labels{end+1,1} = '';
                    end
                end
            end
            paramTable = table(names, values, units, labels, 'VariableNames', {'param','value','unit','label'});
            paramTable = sortrows(paramTable,'param');
            paramStruct = table2struct(paramTable,'ToScalar',true);

            simStruct = struct();
            simStruct.f_value = f_values(i);
            simStruct.results = prdSim;
            simStruct.parameters = paramStruct;
            simStruct.times = times;

            allSim = [allSim; simStruct];          %#ok<AGROW>
        end
    end

    if exist('hScatter', 'var')
        uistack(hScatter,'top');  % bring scatter to front
    end

    xlabel('Time since birth (d)', 'Color', 'k');
    ylabel('Wet weight (g)', 'Color', 'k');
    title(sprintf('%s - %s variation', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color', 'k');

    set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
    grid on

    colormap(cmap)
    cb = colorbar;
    cb.Label.String = strrep(f_mode,'_','\_');
    cb.Color = 'k';
    caxis([f_min f_max])

    lgd = legend('TextColor', 'k', 'Location', 'best');
    set(lgd, 'Color', 'none', 'Box', 'off'); 
    hold off;

    %% 8) Save all simulations + parameters in JSON
    jsonStr = jsonencode(allSim);
    jsonFilename = fullfile(res_dir, sprintf('simulations_CT36EFSA_%s_%s.json', species_name, f_mode));
    fid = fopen(jsonFilename,'w');
    if fid == -1
        error('Unable to create JSON file: %s', jsonFilename);
    end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
    disp(['JSON successfully saved: ', jsonFilename]);

    %% 9) Save extracted outputs to CSV
    results_table = array2table(all_results, 'VariableNames', vars_to_extract);
    results_table.f_value = f_values';  
    results_table = [results_table(:, end), results_table(:, 1:end-1)];

    csv_filename = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_mode));
    writetable(results_table, csv_filename);
    disp(['Outputs saved in : ', csv_filename]);

    %% 10) Save figure as PNG and JPEG
    figFilenamePNG = fullfile(res_dir, sprintf('%s_%s.png', species_name, f_mode));
    figFilenameJPG = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_mode));

    set(gcf,'InvertHardcopy','off');
    saveas(gcf, figFilenamePNG);
    saveas(gcf, figFilenameJPG);

    disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);

    %% 11) Plot extracted variables vs f_value
    figure( ...
        'Name', sprintf('Variables vs f_value (%s)', strrep(f_mode,'_',' ')), ...
        'Color','w', 'Position', [100, 100, 1200, 700]); 

    n_vars = length(vars_to_extract);
    n_cols = ceil(n_vars / 2);

    t = tiledlayout(2, n_cols, 'TileSpacing', 'normal', 'Padding', 'compact');
    xlabel(t, sprintf('Reducing factor applied to %s base value', strrep(f_mode,'_',' ')), ...
        'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

    for v = 1:n_vars
        nexttile;
        plot(f_values, all_results(:, v), '-o', ...
            'Color', [0.2 0.5 0.8], ...
            'MarkerFaceColor', [0.2 0.5 0.8], ...
            'MarkerEdgeColor', [0.2 0.5 0.8], ...
            'LineWidth', 2, ...
            'MarkerSize', 6);

        title(vars_to_extract{v}, 'FontWeight', 'bold', 'FontSize', 14, 'Color', 'k');
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
        grid off; box on;
    end

    figVarsPNG = fullfile(res_dir, sprintf('%s_%s_vars.png', species_name, f_mode));
    figVarsJPG = fullfile(res_dir, sprintf('%s_%s_vars.jpg', species_name, f_mode));
    set(gcf,'InvertHardcopy','off');
    saveas(gcf, figVarsPNG);
    saveas(gcf, figVarsJPG);
    disp(['Variables figure successfully saved: ', figVarsPNG, ' and ', figVarsJPG]);
end

end