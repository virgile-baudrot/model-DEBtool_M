% master_run.m pour Perdix_perdix
% Script d'automatisation pour lancer les prédictions, combiner les images 
% et générer un graphique comparatif des modalités.

clear; close all; clc;

% --- Paramètres ---
species_name = 'Perdix_perdix';
res_dir = sprintf('results_%s', species_name);
f_modes = {'f_nest', 'f_juvenile', 'f_adult'};

%% 1. LANCER LES SIMULATIONS
disp('==================================================');
disp('1. LANCEMENT DES SIMULATIONS DEB POUR PERDIX PERDIX');
disp('==================================================');
% Appelle le script principal de prédiction
run_prediction_Perdix_CT36EFSA(f_modes);


%% 2. COMBINER LES IMAGES DE CROISSANCE
disp('==================================================');
disp('2. COMBINAISON DES GRAPHIQUES DE CROISSANCE');
disp('==================================================');

img_nest = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{1}));
img_juv  = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{2}));
img_ad   = fullfile(res_dir, sprintf('%s_%s.jpg', species_name, f_modes{3}));

if exist(img_nest, 'file') && exist(img_juv, 'file') && exist(img_ad, 'file')
    % Lecture et concaténation horizontale
    img_combinee = [imread(img_nest), imread(img_juv), imread(img_ad)];
    
    % Sauvegarde
    combined_name = fullfile(res_dir, sprintf('%s_croissance_combinee.jpg', species_name));
    imwrite(img_combinee, combined_name);
    disp(['--> Image combinée avec succès : ', combined_name]);
else
    warning('Certaines images de croissance sont introuvables. La concaténation a été ignorée.');
end


%% 3. CREER LE GRAPHIQUE COMPARATIF DES PARAMETRES (FACETTING)
disp('==================================================');
disp('3. CREATION DU GRAPHIQUE COMPARATIF DES VARIABLES');
disp('==================================================');

% Noms des fichiers CSV générés à l'étape 1
csv_nest = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{1}));
csv_juv  = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{2}));
csv_ad   = fullfile(res_dir, sprintf('outputs_%s_%s.csv', species_name, f_modes{3}));

% Lecture des tables
T_nest = readtable(csv_nest);
T_juv  = readtable(csv_juv);
T_ad   = readtable(csv_ad);

% Récupération des noms de variables (en ignorant la colonne 'f_value')
vars_to_extract = T_nest.Properties.VariableNames;
vars_to_extract(strcmp(vars_to_extract, 'f_value')) = []; 

n_vars = length(vars_to_extract);
n_cols = ceil(n_vars / 2);

% --- Dictionnaire des descriptions des variables ---
var_descriptions = containers.Map(...
    {'ab', 'tx', 'tp', 'tR', 'am', 'Wwb', 'Wwi', 'Wwim', 'Ww0', 'Ri'}, ...
    {'age at birth', ...
     'time since birth at fledging', ...
     'time since birth at puberty', ...
     'time since birth at 1st brood', ...
     'life span', ...
     'wet weight at birth', ...
     'ultimate wet weight', ...
     'ultimate wet weight for males', ...
     'initial wet weight', ...
     'maximum reproductive rate'});

% Création de la figure
fig = figure('Name', 'Comparaison des Modalites', 'Color', 'w', 'Position', [100, 100, 1400, 400]);
t = tiledlayout(2, n_cols, 'TileSpacing', 'normal', 'Padding', 'compact');
xlabel(t, 'Reduction factor applied to functional response', 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');

% Définition des couleurs pour les modalités (Orange, Vert, Bleu)
c_nest = [0.8500, 0.3250, 0.0980];
c_juv  = [0.4660, 0.6740, 0.1880];
c_ad   = [0.0000, 0.4470, 0.7410];

for v = 1:n_vars
    var_name = vars_to_extract{v};
    nexttile;
    hold on;

    % Tracer les 3 lignes pour la variable en cours
    p1 = plot(T_nest.f_value, T_nest.(var_name), '-o', 'Color', c_nest, 'MarkerFaceColor', c_nest, 'LineWidth', 2, 'MarkerSize', 4);
    p2 = plot(T_juv.f_value,  T_juv.(var_name),  '-s', 'Color', c_juv,  'MarkerFaceColor', c_juv,  'LineWidth', 2, 'MarkerSize', 4);
    p3 = plot(T_ad.f_value,   T_ad.(var_name),   '-^', 'Color', c_ad,   'MarkerFaceColor', c_ad,   'LineWidth', 2, 'MarkerSize', 4);

    % --- Formatage du titre sur plusieurs lignes ---
    if isKey(var_descriptions, var_name)
        % Création d'un "cell array" pour forcer le retour à la ligne
        % Ligne 1: description complète, Ligne 2: (nom_variable)
        plot_title = {var_descriptions(var_name), sprintf('(%s)', var_name)};
    else
        plot_title = var_name; % Sécurité au cas où une variable manque dans le dictionnaire
    end

    % Title and facet style (FontSize légèrement réduit à 12 pour tenir sur 2 lignes)
    title(plot_title, 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k');
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11);
    grid off; box on;
    hold off;

    % Add legend only to the first plot to avoid cluttering
    if v == 1
        legend([p1, p2, p3], strrep(f_modes, '_', '\_'), 'Location', 'best', 'FontSize', 11);
    end
end

% Sauvegarde de la figure
figCombPNG = fullfile(res_dir, sprintf('%s_comparaison_variables.png', species_name));
figCombJPG = fullfile(res_dir, sprintf('%s_comparaison_variables.jpg', species_name));
set(gcf,'InvertHardcopy','off');
saveas(gcf, figCombPNG);
saveas(gcf, figCombJPG);

disp(['--> Graphique des variables combinées sauvegardé : ', figCombJPG]);
disp('TERMINÉ AVEC SUCCÈS !');