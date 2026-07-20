% master_run.m
% Automation script to run predictions, combine images 
% and generate a comparative plot of modalities.

clear; close all; clc;

% --- Parameters ---
species_name = 'Sorex_araneus';
res_dir = sprintf('results_%s', species_name);
f_modes = {'f_lactation', 'f_juvenile', 'f_adult'};

%% 1. RUN SIMULATIONS
disp('==================================================');
disp('1. RUNNING DEB SIMULATIONS');
disp('==================================================');
run_prediction_Sorex_araneus_CT36EFSA(f_modes);


%% 2. COMBINE GROWTH IMAGES
disp('==================================================');
disp('2. COMBINING GROWTH PLOTS');
disp('==================================================');

img_lact = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{1}));
img_juv  = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{2}));
img_ad   = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{3}));

if exist(img_lact, 'file') && exist(img_juv, 'file') && exist(img_ad, 'file')
    % Read and concatenate horizontally
    img_combined = [imread(img_lact), imread(img_juv), imread(img_ad)];
    % Save
    combined_name = fullfile(res_dir, sprintf('%s_combined_growth.jpg', species_name));
    imwrite(img_combined, combined_name);
    disp(['--> Combined image saved: ', combined_name]);
else
    warning('Some growth images are missing. Concatenation skipped.');
end


%% 3. CREATE COMPARATIVE PLOT (FACETTING)
disp('==================================================');
disp('3. CREATING COMPARATIVE VARIABLES PLOT');
disp('==================================================');

% CSV files generated in step 1
csv_lact = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{1}));
csv_juv  = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{2}));
csv_ad   = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{3}));

% Read tables
T_lact = readtable(csv_lact);
T_juv  = readtable(csv_juv);
T_ad   = readtable(csv_ad);

% Get variable names (ignoring 'f_value' column)
vars_to_extract = T_lact.Properties.VariableNames;
vars_to_extract(strcmp(vars_to_extract, 'f_value')) = []; 

n_vars = length(vars_to_extract);
n_cols = ceil(n_vars / 2);

% --- Map of variables descriptions (Mammals) ---
var_descriptions = containers.Map(...
    {'tg', 'tx', 'tp', 'am', 'Li', 'Wwb', 'Wwx', 'Wwp', 'Wwi', 'Ri'}, ...
    {'gestation time', ...
     'time since birth at weaning', ...
     'time since birth at puberty', ...
     'life span', ...
     'ultimate total length', ...
     'wet weight at birth', ...
     'wet weight at weaning', ...
     'wet weight at puberty', ...
     'ultimate wet weight', ...
     'maximum reproductive rate'});

% Create figure (Height reduced to 400)
fig = figure('Name', 'Modalities Comparison', 'Color', 'w', 'Position', [100, 100, 1400, 400]);
t = tiledlayout(2, n_cols, 'TileSpacing', 'normal', 'Padding', 'compact');
xlabel(t, 'Factor applied to function response value', 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

% Define colors for modalities (Orange, Green, Blue)
c_lact = [0.8500, 0.3250, 0.0980];
c_juv  = [0.4660, 0.6740, 0.1880];
c_ad   = [0.0000, 0.4470, 0.7410];

for v = 1:n_vars
    var_name = vars_to_extract{v};
    nexttile;
    hold on;

    % Plot 3 lines for current variable
    p1 = plot(T_lact.f_value, T_lact.(var_name), '-o', 'Color', c_lact, 'MarkerFaceColor', c_lact, 'LineWidth', 2, 'MarkerSize', 4);
    p2 = plot(T_juv.f_value,  T_juv.(var_name),  '-s', 'Color', c_juv,  'MarkerFaceColor', c_juv,  'LineWidth', 2, 'MarkerSize', 4);
    p3 = plot(T_ad.f_value,   T_ad.(var_name),   '-^', 'Color', c_ad,   'MarkerFaceColor', c_ad,   'LineWidth', 2, 'MarkerSize', 4);

    % --- Multi-line Title Formatting ---
    if isKey(var_descriptions, var_name)
        % Create a 2-line title using a cell array
        plot_title = {var_descriptions(var_name), sprintf('(%s)', var_name)};
    else
        plot_title = var_name; % Fallback if variable is not in the Map
    end

    % Title and facet style
    title(plot_title, 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11);
    grid off; box on;
    hold off;

    % Add legend only to the first plot to avoid cluttering
    if v == 1
        legend([p1, p2, p3], strrep(f_modes, '_', '\_'), 'Location', 'best', 'FontSize', 11);
    end
end

% Save figure
figCombPNG = fullfile(res_dir, sprintf('%s_variables_comparison.png', species_name));
figCombJPG = fullfile(res_dir, sprintf('%s_variables_comparison.jpg', species_name));
set(gcf,'InvertHardcopy','off');
saveas(gcf, figCombPNG);
saveas(gcf, figCombJPG);

disp(['--> Combined variables plot saved: ', figCombJPG]);
disp('FINISHED SUCCESSFULLY!');