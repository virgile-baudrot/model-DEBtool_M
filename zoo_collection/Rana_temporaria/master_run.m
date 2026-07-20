% master_run.m
% Automation script to run predictions, combine images 
% and generate a comparative plot of modalities.

clear; close all; clc;

% --- Parameters ---
species_name = 'Rana_temporaria';
res_dir = sprintf('results_%s', species_name);
f_modes = {'f_nest', 'f_juvenile', 'f_adult'};

%% 1. RUN SIMULATIONS
disp('==================================================');
disp('1. RUNNING DEB SIMULATIONS');
disp('==================================================');
run_prediction_Rana_temporaria_CT36EFSA(f_modes);


%% 2. COMBINE GROWTH IMAGES (FEMALES, MALES, LdL)
disp('==================================================');
disp('2. COMBINING GROWTH PLOTS');
disp('==================================================');

% Load Female plots
img_nest_f = fullfile(res_dir, sprintf('%s_%s_female.jpg', species_name, f_modes{1}));
img_juv_f  = fullfile(res_dir, sprintf('%s_%s_female.jpg', species_name, f_modes{2}));
img_ad_f   = fullfile(res_dir, sprintf('%s_%s_female.jpg', species_name, f_modes{3}));

% Load Male plots
img_nest_m = fullfile(res_dir, sprintf('%s_%s_male.jpg', species_name, f_modes{1}));
img_juv_m  = fullfile(res_dir, sprintf('%s_%s_male.jpg', species_name, f_modes{2}));
img_ad_m   = fullfile(res_dir, sprintf('%s_%s_male.jpg', species_name, f_modes{3}));

% Load Change in Length (LdL) plots
img_nest_LdL = fullfile(res_dir, sprintf('%s_%s_LdL.jpg', species_name, f_modes{1}));
img_juv_LdL  = fullfile(res_dir, sprintf('%s_%s_LdL.jpg', species_name, f_modes{2}));
img_ad_LdL   = fullfile(res_dir, sprintf('%s_%s_LdL.jpg', species_name, f_modes{3}));

if exist(img_nest_f, 'file') && exist(img_ad_m, 'file') && exist(img_ad_LdL, 'file')
    % Concatenate SVL (Time vs Length) for Females and Males (2 rows x 3 cols)
    img_combined_time = [imread(img_nest_f), imread(img_juv_f), imread(img_ad_f) ; ...
                         imread(img_nest_m), imread(img_juv_m), imread(img_ad_m)];
    
    % Concatenate Change in Length (1 row x 3 cols)
    img_combined_LdL = [imread(img_nest_LdL), imread(img_juv_LdL), imread(img_ad_LdL)];
    
    % Save Both
    combined_time_name = fullfile(res_dir, sprintf('%s_combined_growth_time.jpg', species_name));
    combined_LdL_name = fullfile(res_dir, sprintf('%s_combined_growth_LdL.jpg', species_name));
    
    imwrite(img_combined_time, combined_time_name);
    imwrite(img_combined_LdL, combined_LdL_name);
    disp(['--> Combined Time Growth image saved: ', combined_time_name]);
    disp(['--> Combined LdL image saved: ', combined_LdL_name]);
else
    warning('Some growth images are missing. Concatenation skipped.');
end


%% 3. CREATE COMPARATIVE PLOT (FACETTING)
disp('==================================================');
disp('3. CREATING COMPARATIVE VARIABLES PLOT');
disp('==================================================');

% CSV files generated in step 1
csv_nest = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{1}));
csv_juv  = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{2}));
csv_ad   = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{3}));

% Read tables
T_nest = readtable(csv_nest);
T_juv  = readtable(csv_juv);
T_ad   = readtable(csv_ad);

% Get variable names (ignoring 'f_value' column)
vars_to_extract = T_nest.Properties.VariableNames;
vars_to_extract(strcmp(vars_to_extract, 'f_value')) = []; 

n_vars = length(vars_to_extract);
n_cols = ceil(n_vars / 2);

% --- Map of variables descriptions (Amphibians) ---
var_descriptions = containers.Map(...
    {'am', 'tj', 'tp', 'Lj', 'Lp', 'Lpm', 'Li', 'Lim', 'Wwb', 'Wwi', 'Wwim', 'Ri'}, ...
    {'life span', ...
     'time since birth at metamorphosis', ...
     'time since birth at puberty', ...
     'SVL at metamorphosis', ...
     'SVL at puberty (females)', ...
     'SVL at puberty (males)', ...
     'ultimate SVL (females)', ...
     'ultimate SVL (males)', ...
     'wet weight at birth', ...
     'ultimate wet weight (females)', ...
     'ultimate wet weight (males)', ...
     'maximum reproductive rate'});

% Create figure (Adjusted height and width for 12 variables)
fig = figure('Name', 'Modalities Comparison', 'Color', 'w', 'Position', [100, 100, 1600, 500]);
t = tiledlayout(2, n_cols, 'TileSpacing', 'normal', 'Padding', 'compact');
xlabel(t, 'Factor applied to function response value', 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

% Define colors for modalities (Orange, Green, Blue)
c_nest = [0.8500, 0.3250, 0.0980];
c_juv  = [0.4660, 0.6740, 0.1880];
c_ad   = [0.0000, 0.4470, 0.7410];

for v = 1:n_vars
    var_name = vars_to_extract{v};
    nexttile;
    hold on;

    % Plot 3 lines for current variable
    p1 = plot(T_nest.f_value, T_nest.(var_name), '-o', 'Color', c_nest, 'MarkerFaceColor', c_nest, 'LineWidth', 2, 'MarkerSize', 4);
    p2 = plot(T_juv.f_value,  T_juv.(var_name),  '-s', 'Color', c_juv,  'MarkerFaceColor', c_juv,  'LineWidth', 2, 'MarkerSize', 4);
    p3 = plot(T_ad.f_value,   T_ad.(var_name),   '-^', 'Color', c_ad,   'MarkerFaceColor', c_ad,   'LineWidth', 2, 'MarkerSize', 4);

    % --- Multi-line Title Formatting ---
    if isKey(var_descriptions, var_name)
        plot_title = {var_descriptions(var_name), sprintf('(%s)', var_name)};
    else
        plot_title = var_name; 
    end

    % Title and facet style
    title(plot_title, 'FontWeight', 'bold', 'FontSize', 11, 'Color', 'k');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11);
    grid off; box on;
    hold off;

    % Add legend only to the first plot to avoid cluttering
    if v == 1
        legend([p1, p2, p3], strrep(f_modes, '_', '\_'), 'Location', 'best', 'FontSize', 10);
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