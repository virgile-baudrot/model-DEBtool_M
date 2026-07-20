% master_run.m
% Automation script to run predictions, combine images 
% and generate a comparative plot of modalities.

clear; close all; clc;

% --- Parameters ---
species_name = 'Podarcis_muralis';
res_dir = sprintf('results_%s', species_name);
f_modes = {'f_embryo', 'f_juvenile', 'f_adult'};

%% 1. RUN SIMULATIONS
disp('==================================================');
disp('1. RUNNING DEB SIMULATIONS');
disp('==================================================');
run_prediction_Podarcis_muralis_CT36EFSA(f_modes);


%% 2. COMBINE GROWTH IMAGES (DYNAMICS & FECUNDITY)
disp('==================================================');
disp('2. COMBINING GROWTH PLOTS');
disp('==================================================');

% --- 2A. COMBINE DYNAMICS: 4 plots/row (Female tL, Male tL, Female LW, Male LW) x 2 rows (Juv, Adult)
dyn_files = { ...
    'f_juvenile_tL_female', 'f_juvenile_tL_male', 'f_juvenile_LW_female', 'f_juvenile_LW_male'; ...
    'f_adult_tL_female',    'f_adult_tL_male',    'f_adult_LW_female',    'f_adult_LW_male' ...
};

dyn_imgs = cell(2, 4);
all_dyn_exist = true;
ref_size = [];

for r = 1:2
    for c = 1:4
        pth = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, dyn_files{r,c}));
        if exist(pth, 'file')
            [im, map] = imread(pth);
            if ~isempty(map), im = im2uint8(ind2rgb(im, map)); end
            if size(im, 3) == 1, im = repmat(im, [1 1 3]); end
            
            if isempty(ref_size)
                ref_size = [size(im, 1), size(im, 2)];
            else
                im = imresize(im, ref_size);
            end
            dyn_imgs{r,c} = im;
        else
            all_dyn_exist = false;
        end
    end
end

if all_dyn_exist
    img_combined_dyn = [dyn_imgs{1,1}, dyn_imgs{1,2}, dyn_imgs{1,3}, dyn_imgs{1,4} ; ...
                        dyn_imgs{2,1}, dyn_imgs{2,2}, dyn_imgs{2,3}, dyn_imgs{2,4}];
    combined_dyn_name = fullfile(res_dir, sprintf('%s_combined_dynamics.jpg', species_name));
    imwrite(img_combined_dyn, combined_dyn_name);
    disp(['--> Combined dynamics image saved: ', combined_dyn_name]);
else
    warning('Some dynamics images are missing. Concatenation skipped.');
end

% --- 2B. COMBINE FECUNDITY: 3 plots/row (Emb, Juv, Ad) x 1 row
fec_files = {'f_embryo_LN', 'f_juvenile_LN', 'f_adult_LN'};
fec_imgs = cell(1, 3);
all_fec_exist = true;

for c = 1:3
    pth = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, fec_files{c}));
    if exist(pth, 'file')
        [im, map] = imread(pth);
        if ~isempty(map), im = im2uint8(ind2rgb(im, map)); end
        if size(im, 3) == 1, im = repmat(im, [1 1 3]); end
        
        if ~isempty(ref_size)
            im = imresize(im, ref_size);
        end
        fec_imgs{c} = im;
    else
        all_fec_exist = false;
    end
end

if all_fec_exist
    img_combined_fec = [fec_imgs{1}, fec_imgs{2}, fec_imgs{3}];
    combined_fec_name = fullfile(res_dir, sprintf('%s_combined_fecundity.jpg', species_name));
    imwrite(img_combined_fec, combined_fec_name);
    disp(['--> Combined fecundity image saved: ', combined_fec_name]);
else
    warning('Some fecundity images are missing. Concatenation skipped.');
end


%% 3. CREATE COMPARATIVE PLOT (FACETTING)
disp('==================================================');
disp('3. CREATING COMPARATIVE VARIABLES PLOT');
disp('==================================================');

% CSV files generated in step 1
csv_emb = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{1}));
csv_juv = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{2}));
csv_ad  = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{3}));

% Read tables
T_emb = readtable(csv_emb);
T_juv = readtable(csv_juv);
T_ad  = readtable(csv_ad);

% Get variable names (ignoring 'f_value' column)
vars_to_extract = T_emb.Properties.VariableNames;
vars_to_extract(strcmp(vars_to_extract, 'f_value')) = []; 

n_vars = length(vars_to_extract);
n_cols = ceil(n_vars / 2);

% --- Map of variables descriptions (Reptiles) ---
var_descriptions = containers.Map(...
    {'tp', 'am', 'Lb', 'Lp', 'Lpm', 'Li', 'Lim', 'Ww0', 'Wwb', 'Wwi', 'Wwim', 'Ri'}, ...
    {'time since birth at puberty', ...
     'life span', ...
     'SVL at birth', ...
     'SVL at puberty (females)', ...
     'SVL at puberty (males)', ...
     'ultimate SVL (females)', ...
     'ultimate SVL (males)', ...
     'wet weight of eggs', ...
     'wet weight at birth', ...
     'ultimate wet weight (females)', ...
     'ultimate wet weight (males)', ...
     'maximum reproductive rate'});

% Create figure 
fig = figure('Name', 'Modalities Comparison', 'Color', 'w', 'Position', [100, 100, 1600, 500]);
t = tiledlayout(2, n_cols, 'TileSpacing', 'normal', 'Padding', 'compact');
xlabel(t, 'Factor applied to function response value', 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

% Define colors for modalities (Orange, Green, Blue)
c_emb = [0.8500, 0.3250, 0.0980];
c_juv = [0.4660, 0.6740, 0.1880];
c_ad  = [0.0000, 0.4470, 0.7410];

for v = 1:n_vars
    var_name = vars_to_extract{v};
    nexttile;
    hold on;

    % Plot 3 lines for current variable
    p1 = plot(T_emb.f_value, T_emb.(var_name), '-o', 'Color', c_emb, 'MarkerFaceColor', c_emb, 'LineWidth', 2, 'MarkerSize', 4);
    p2 = plot(T_juv.f_value, T_juv.(var_name), '-s', 'Color', c_juv, 'MarkerFaceColor', c_juv, 'LineWidth', 2, 'MarkerSize', 4);
    p3 = plot(T_ad.f_value,  T_ad.(var_name),  '-^', 'Color', c_ad,  'MarkerFaceColor', c_ad,  'LineWidth', 2, 'MarkerSize', 4);

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

% Save figure robustly
figCombPNG = fullfile(res_dir, sprintf('%s_variables_comparison.png', species_name));
figCombJPG = fullfile(res_dir, sprintf('%s_variables_comparison.jpg', species_name));
exportgraphics(fig, figCombPNG, 'BackgroundColor', 'white', 'Resolution', 150);
exportgraphics(fig, figCombJPG, 'BackgroundColor', 'white', 'Resolution', 150);

disp(['--> Combined variables plot saved: ', figCombJPG]);
disp('FINISHED SUCCESSFULLY!');