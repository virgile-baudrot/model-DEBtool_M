%% === Perdix perdix : prédictions DEB + données (points blancs) + tableau des paramètres ===

clear; close all; clc;

%% 1) Charger données & paramètres
[data, auxData, metaData, txtData, weights] = mydata_Perdix_perdix();
[par, metaPar, txtPar] = pars_init_Perdix_perdix_discrete(metaData);

%% 2) Préparer le vecteur de temps pour la simulation
t_max  = 150;                          % durée de simulation en jours
times  = linspace(0, t_max, 301)';      % 301 points réguliers
dataSim = data;                         
dataSim.tW = [times, zeros(size(times))]; % 2e colonne placeholder

%% 3) Choisir les valeurs de f à comparer
f_values = [1.0, 1.1, 0.9];  % ad libitum, modéré, limité

%% 4) Tracer les prédictions + points expérimentaux
figure('Name','Growth of Perdix perdix','Color','k'); hold on; box on;
colors = lines(numel(f_values));

% --- Données expérimentales : points blancs avec contour noir
if isfield(data,'tW') && ~isempty(data.tW)
    scatter(data.tW(:,1), data.tW(:,2), 50, 'w', 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',1.2, ...
        'DisplayName','Experimental data');
end

% --- Boucle sur les différentes valeurs de f
for i = 1:numel(f_values)
    par.f_juvenile = f_values(i);
    [prdSim, info] = predict_Perdix_perdix_discrete(par, dataSim, auxData);
    if info ~= 0 && isfield(prdSim,'tW')
        plot(times, prdSim.tW, 'LineWidth', 2, ...
            'Color', colors(i,:), ...
            'DisplayName', sprintf('Prediction (f = %.1f)', f_values(i)));
    end
end

% --- Style graphique
xlabel('time since brith (d)','Color','w');
ylabel('wet weight (g)','Color','w');
title('Perdix perdix - Growth (predictions vs data)','Color','w');
legend('TextColor','w','Location','best');
set(gca,'Color','k','XColor','w','YColor','w','FontSize',11);
grid on; hold off;

%% 5) Construire le tableau des paramètres (valeur + unité + libellé)
allNames = fieldnames(par);
names = {}; values = []; units = {}; labels = {};

for k = 1:numel(allNames)
    nm = allNames{k};
    if strcmp(nm,'free'); continue; end
    val = par.(nm);
    if isnumeric(val) && isscalar(val)
        names{end+1,1} = nm; %#ok<SAGROW>
        values(end+1,1) = val; %#ok<SAGROW>
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

paramTable = table(names, values, units, labels, ...
    'VariableNames', {'param','value','unit','label'});
paramTable = sortrows(paramTable,'param');
