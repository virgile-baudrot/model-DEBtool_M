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

paths_time = {img_nest_f, img_juv_f, img_ad_f, img_nest_m, img_juv_m, img_ad_m};
paths_LdL  = {img_nest_LdL, img_juv_LdL, img_ad_LdL};

% --- ROBUST IMAGE READING FOR TIME PLOTS ---
imgs_time = cell(1, 6);
all_time_exist = true;
for k = 1:6
    if exist(paths_time{k}, 'file')
        [im, map] = imread(paths_time{k});
        if ~isempty(map)
            im = im2uint8(ind2rgb(im, map));
        end
        if size(im, 3) == 1
            im = repmat(im, [1 1 3]); 
        end
        imgs_time{k} = im;
    else
        all_time_exist = false;
    end
end

% --- ROBUST IMAGE READING FOR LDL PLOTS ---
imgs_LdL = cell(1, 3);
all_LdL_exist = true;
for k = 1:3
    if exist(paths_LdL{k}, 'file')
        [im, map] = imread(paths_LdL{k});
        if ~isempty(map)
            im = im2uint8(ind2rgb(im, map));
        end
        if size(im, 3) == 1
            im = repmat(im, [1 1 3]); 
        end
        imgs_LdL{k} = im;
    else
        all_LdL_exist = false;
    end
end

if all_time_exist
    ref_size_time = [size(imgs_time{1}, 1), size(imgs_time{1}, 2)];
    for k = 2:6
        imgs_time{k} = imresize(imgs_time{k}, ref_size_time);
    end
    img_combined_time = [imgs_time{1}, imgs_time{2}, imgs_time{3} ; ...
                         imgs_time{4}, imgs_time{5}, imgs_time{6}];
    
    combined_time_name = fullfile(res_dir, sprintf('%s_combined_growth_time.jpg', species_name));
    imwrite(img_combined_time, combined_time_name);
    disp(['--> Combined Time Growth image saved: ', combined_time_name]);
else
    warning('Some Time growth images are missing. Concatenation skipped.');
end

if all_LdL_exist
    ref_size_LdL = [size(imgs_LdL{1}, 1), size(imgs_LdL{1}, 2)];
    for k = 2:3
        imgs_LdL{k} = imresize(imgs_LdL{k}, ref_size_LdL);
    end
    img_combined_LdL = [imgs_LdL{1}, imgs_LdL{2}, imgs_LdL{3}];
    
    combined_LdL_name = fullfile(res_dir, sprintf('%s_combined_growth_LdL.jpg', species_name));
    imwrite(img_combined_LdL, combined_LdL_name);
    disp(['--> Combined LdL image saved: ', combined_LdL_name]);
else
    warning('Some LdL images are missing. Concatenation skipped.');
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

% Robust Export
figCombPNG = fullfile(res_dir, sprintf('%s_variables_comparison.png', species_name));
figCombJPG = fullfile(res_dir, sprintf('%s_variables_comparison.jpg', species_name));
exportgraphics(fig, figCombPNG, 'BackgroundColor', 'white', 'Resolution', 150);
exportgraphics(fig, figCombJPG, 'BackgroundColor', 'white', 'Resolution', 150);

disp(['--> Combined variables plot saved: ', figCombJPG]);
disp('FINISHED SUCCESSFULLY!');