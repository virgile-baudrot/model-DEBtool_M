%% === DEB predictions + data + save JSON simulations and parameters ===
clear; close all; clc;

%% 1) Species name (can be changed to any other species)
%% 4) Choose the f values to simulate
species_name = "Podarcis_muralis"; t_max  = 6000; f_values = 0.5:0.01:1.1;

% --- Which f is varied in the simulations
% Possible values: "f_juvenile", "f_adult"
f_mode = "f_juvenile";
% f_mode = "f_adult";

add_shift = -0.1;  % is initial f_xxx -1 - default is 0

%% --- Generate dynamic function names based on species_name
dataFunc    = str2func(sprintf('mydata_%s', species_name));
parsFunc    = str2func(sprintf('pars_init_%s_CT36EFSA', species_name));
predictFunc = str2func(sprintf('predict_%s_CT36EFSA', species_name));

%% 2) Load data and parameters
[data, auxData, metaData, txtData, weights] = dataFunc();
[par, metaPar, txtPar] = parsFunc(metaData);

%% 3) Prepare the time vector for simulation
%                         % simulation duration in days
times  = linspace(0, t_max, 301)';      % 301 equally spaced points
SVLs  = linspace(5, 7, 301)';      % 301 equally spaced points
dataSim = data;
dataSim.tL = [times, zeros(size(times))]; % placeholder for simulated weights
dataSim.SVLs = [SVLs, zeros(size(SVLs))]; % placeholder for simulated weights

%% AFTER THAT LINE, EVERYTHING IS GENERIC

%% 5) Define color mapping (like viridis/plasma)
f_min = 0.5; f_max = 1.1; df = 0.01;
f_ref = f_min:df:f_max; nCol = numel(f_ref);
cmap = flipud(turbo(nCol));             % flip turbo so blue = high f
f2idx = @(f) max(1, min(nCol, round((f - f_min)/df) + 1)); % map f -> color index

%% 6) Plot predictions + experimental data
% --- Loop over f values and plot predictions
alpha_val = 0.7;      % transparency for non-max curves
lw_default = 2;       % default line width
lw_f1 = 3;            % line width for f = 1
%% GRAPHIC 1 - FEMALE
figure( ...
    'Name', sprintf('SVL of %s female', strrep(species_name,'_',' ')), ...
    'Color','w');   % fond blanc à l'écran
hold on; box on;
% --- Scatter experimental data on top
if isfield(data,'tL_D') && ~isempty(data.tL_D)
    hScatter = scatter(data.tL_D(:,1), data.tL_D(:,2), 50, 'w', 'filled', ...
                   'MarkerEdgeColor','k', 'LineWidth',1.2);
end
for i = 1:numel(f_values)
    par.(f_mode) = f_values(i) + add_shift;
    [prdSim, info] = predictFunc(par, dataSim, auxData);

    if info ~= 0 && isfield(prdSim,'tL_D')
        idx = f2idx(f_values(i));
        % idx = f2idx(par.(f_mode));

        if abs(par.(f_mode) - 1.0) < 1e-6
            % --- f = 1: opaque + wider line
            plot(times, prdSim.tL_D(:,1), 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
        else
            % --- other f: semi-transparent
            plot(times, prdSim.tL_D(:,1), 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
        end
    end
end
uistack(hScatter,'top');  % force le scatter en premier plan
% --- Axis labels and title
xlabel('time since birth (d)','Color','k');
ylabel('SVL (cm)','Color','k');
title(sprintf('%s female - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color','k');

% --- Customize axes
set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
grid on

% --- Colormap and colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = strrep(f_mode,'_','\_');
cb.Color = 'k';
caxis([f_min f_max])

legend('TextColor','k','Location','best'); 
legend off
hold off;

%% --- Save figure as PNG and JPEG with species name
figFilenamePNG = sprintf('%s_%s_female.png', species_name, f_mode);
figFilenameJPG = sprintf('%s_%s_female.jpg', species_name, f_mode);

set(gcf,'InvertHardcopy','off');
saveas(gcf, figFilenamePNG);
saveas(gcf, figFilenameJPG);

disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);

%% GRAPHIC 2 - MALE
figure( ...
    'Name', sprintf('SVL of %s male', strrep(species_name,'_',' ')), ...
    'Color','w');   % fond blanc à l'écran
hold on; box on;
% --- Scatter experimental data on top
if isfield(data,'tL_D') && ~isempty(data.tL_D)
    hScatter = scatter(data.tL_D(:,1), data.tL_D(:,3), 50, 'w', 'filled', ...
                   'MarkerEdgeColor','k', 'LineWidth',1.2);
end
for i = 1:numel(f_values)
    par.(f_mode) = f_values(i) + add_shift;
    [prdSim, info] = predictFunc(par, dataSim, auxData);

    if info ~= 0 && isfield(prdSim,'tL_D')
        idx = f2idx(f_values(i));
        % idx = f2idx(par.(f_mode));

        if abs(par.(f_mode) - 1.0) < 1e-6
            % --- f = 1: opaque + wider line
            plot(times, prdSim.tL_D(:,2), 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
        else
            % --- other f: semi-transparent
            plot(times, prdSim.tL_D(:,2), 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
        end
    end
end
uistack(hScatter,'top');  % force le scatter en premier plan
% --- Axis labels and title
xlabel('time since birth (d)','Color','k');
ylabel('SVL (cm)','Color','k');
title(sprintf('%s male - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color','k');

% --- Customize axes
set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
grid on

% --- Colormap and colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = strrep(f_mode,'_','\_');
cb.Color = 'k';
caxis([f_min f_max])

legend('TextColor','k','Location','best'); 
legend off
hold off;


%% --- Save figure as PNG and JPEG with species name
figFilenamePNG = sprintf('%s_%s_male.png', species_name, f_mode);
figFilenameJPG = sprintf('%s_%s_male.jpg', species_name, f_mode);

set(gcf,'InvertHardcopy','off');
saveas(gcf, figFilenamePNG);
saveas(gcf, figFilenameJPG);

disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);

%% GRAPHIC 3 - FEMALE Ww ~ SVL
figure( ...
    'Name', sprintf('Wet weight of %s female', strrep(species_name,'_',' ')), ...
    'Color','w');   % fond blanc à l'écran
hold on; box on;
% --- Scatter experimental data on top
if isfield(data,'LW_f') && ~isempty(data.LW_f)
    hScatter = scatter(data.LW_f(:,1), data.LW_f(:,2), 50, 'w', 'filled', ...
                   'MarkerEdgeColor','k', 'LineWidth',1.2);
end
for i = 1:numel(f_values)
    par.(f_mode) = f_values(i) + add_shift;
    [prdSim, info] = predictFunc(par, dataSim, auxData);

    if info ~= 0 && isfield(prdSim,'LW_f')
        idx = f2idx(f_values(i));
        % idx = f2idx(par.(f_mode));

        if abs(par.(f_mode) - 1.0) < 1e-6
            % --- f = 1: opaque + wider line
            plot(prdSim.tL_D(:,1), prdSim.LW_f, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
        else
            % --- other f: semi-transparent
            plot(prdSim.tL_D(:,1), prdSim.LW_f, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
        end
    end
end
uistack(hScatter,'top');  % force le scatter en premier plan
% --- Axis labels and title
xlabel('SVL (cm)','Color','k');
ylabel('wet weight (g)','Color','k');
title(sprintf('%s female - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color','k');

% --- Customize axes
set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
grid on

% --- Colormap and colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = strrep(f_mode,'_','\_');
cb.Color = 'k';
caxis([f_min f_max])

legend('TextColor','k','Location','best'); 
legend off
hold off;


%% --- Save figure as PNG and JPEG with species name
figFilenamePNG = sprintf('%s_%s_female_ww.png', species_name, f_mode);
figFilenameJPG = sprintf('%s_%s_female_ww.jpg', species_name, f_mode);

set(gcf,'InvertHardcopy','off');
saveas(gcf, figFilenamePNG);
saveas(gcf, figFilenameJPG);

disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);


%% GRAPHIC 4 - MALE Ww ~ SVL
figure( ...
    'Name', sprintf('Wet weight of %s male', strrep(species_name,'_',' ')), ...
    'Color','w');   % fond blanc à l'écran
hold on; box on;
% --- Scatter experimental data on top
if isfield(data,'LW_m') && ~isempty(data.LW_m)
    hScatter = scatter(data.LW_m(:,1), data.LW_m(:,2), 50, 'w', 'filled', ...
                   'MarkerEdgeColor','k', 'LineWidth',1.2);
end
for i = 1:numel(f_values)
    par.(f_mode) = f_values(i) + add_shift;
    [prdSim, info] = predictFunc(par, dataSim, auxData);

    if info ~= 0 && isfield(prdSim,'LW_m')
        idx = f2idx(f_values(i));
        % idx = f2idx(par.(f_mode));

        if abs(par.(f_mode) - 1.0) < 1e-6
            % --- f = 1: opaque + wider line
            plot(prdSim.tL_D(:,2), prdSim.LW_m, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
        else
            % --- other f: semi-transparent
            plot(prdSim.tL_D(:,2), prdSim.LW_m, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
        end
    end
end
uistack(hScatter,'top');  % force le scatter en premier plan
% --- Axis labels and title
xlabel('SVL (cm)','Color','k');
ylabel('wet weight (g)','Color','k');
title(sprintf('%s male - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color','k');

% --- Customize axes
set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
grid on

% --- Colormap and colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = strrep(f_mode,'_','\_');
cb.Color = 'k';
caxis([f_min f_max])

legend('TextColor','k','Location','best'); 
legend off
hold off;


%% --- Save figure as PNG and JPEG with species name
figFilenamePNG = sprintf('%s_%s_male_ww.png', species_name, f_mode);
figFilenameJPG = sprintf('%s_%s_male_ww.jpg', species_name, f_mode);

set(gcf,'InvertHardcopy','off');
saveas(gcf, figFilenamePNG);
saveas(gcf, figFilenameJPG);

disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);

%% GRAPHIC 5 - Yearly fecundity ~ SVL
figure( ...
    'Name', sprintf('Yearly Fecundity of %s', strrep(species_name,'_',' ')), ...
    'Color','w');   % fond blanc à l'écran
hold on; box on;
% --- Scatter experimental data on top
if isfield(data,'LN') && ~isempty(data.LN)
    hScatter = scatter(data.LN(:,1), data.LN(:,2), 50, 'w', 'filled', ...
                   'MarkerEdgeColor','k', 'LineWidth',1.2);
end
for i = 1:numel(f_values)
    par.(f_mode) = f_values(i) + add_shift;
    [prdSim, info] = predictFunc(par, dataSim, auxData);

    if info ~= 0 && isfield(prdSim,'LN')
        idx = f2idx(f_values(i));
        % idx = f2idx(par.(f_mode));

        if abs(par.(f_mode) - 1.0) < 1e-6
            % --- f = 1: opaque + wider line
            plot(prdSim.tL_D(:,1), prdSim.LN, 'LineWidth', lw_f1, 'Color', cmap(idx,:), 'HandleVisibility','off');
        else
            % --- other f: semi-transparent
            plot(prdSim.tL_D(:,1), prdSim.LN, 'LineWidth', lw_default, 'Color', [cmap(idx,:) alpha_val], 'HandleVisibility','off');
        end
    end
end
uistack(hScatter,'top');  % force le scatter en premier plan
% --- Axis labels and title
xlabel('SVL (cm)','Color','k');
ylabel('yearly fecundity (#)','Color','k');
title(sprintf('%s - %s', strrep(species_name,'_',' '), strrep(f_mode,'_',' ')), 'Color','k');

% --- Customize axes
set(gca,'Color','w','XColor','k','YColor','k','FontSize',11);
grid on

% --- Colormap and colorbar
colormap(cmap)
cb = colorbar;
cb.Label.String = strrep(f_mode,'_','\_');
cb.Color = 'k';
caxis([f_min f_max])

legend('TextColor','k','Location','best'); 
legend off
hold off;


%% --- Save figure as PNG and JPEG with species name
figFilenamePNG = sprintf('%s_%s_LN.png', species_name, f_mode);
figFilenameJPG = sprintf('%s_%s_LN.jpg', species_name, f_mode);

set(gcf,'InvertHardcopy','off');
saveas(gcf, figFilenamePNG);
saveas(gcf, figFilenameJPG);

disp(['Figure successfully saved: ', figFilenamePNG, ' and ', figFilenameJPG]);

%% 7) Save all simulations + parameters in JSON (struct array)
allSim = [];  % initialize struct array

for iF = 1:numel(f_values)
    par.(f_mode) = f_values(iF);
    [prdSim, info] = predictFunc(par, dataSim, auxData);
    
    if info ~= 0
        %% --- Construct parameter table for this simulation
        allNames = fieldnames(par);
        names = {}; values = []; units = {}; labels = {};
        for k = 1:numel(allNames)
            nm = allNames{k};
            if strcmp(nm,'free'); continue; end
            val = par.(nm);
            if isnumeric(val) && isscalar(val)
                names{end+1,1} = nm; 
                values(end+1,1) = val;
                % --- units
                if isfield(txtPar,'units') && isfield(txtPar.units,nm)
                    units{end+1,1} = txtPar.units.(nm);
                else
                    units{end+1,1} = '';
                end
                % --- labels
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
        
        %% --- Add this simulation to the struct array
        simStruct = struct();
        simStruct.f_value = f_values(iF);      % f_juvenile value
        simStruct.results = prdSim;            % simulation results
        simStruct.parameters = paramStruct;    % associated parameters
        
        allSim = [allSim; simStruct];          %#ok<AGROW>
    end
end


% --- Encode the struct array as JSON
jsonStr = jsonencode(allSim);

%% --- Save struct array to JSON with species name in the filename
jsonFilename = sprintf('simulations_CT36EFSA_%s_%s_.json', species_name, f_mode);

fid = fopen(jsonFilename,'w');
if fid == -1
    error('Unable to create JSON file: %s', jsonFilename);
end
fwrite(fid, jsonStr, 'char'); 
fclose(fid);

disp(['JSON successfully saved: ', jsonFilename]);

