%% Plot fitted roughness from All_40_Dimpled_Corrected_Results.xlsx
% Place this script anywhere; workbook is expected in your Downloads folder.
clear; clc; close all;

filePath = fullfile(getenv('USERPROFILE'),'Downloads', ...
    'All_40_Dimpled_Corrected_Results.xlsx');
% macOS/Linux fallback:
if ~isfile(filePath)
    filePath = fullfile(getenv('HOME'),'Downloads', ...
        'All_40_Dimpled_Corrected_Results.xlsx');
end
if ~isfile(filePath)
    [fileName,folder] = uigetfile('*.xlsx','Select the 40-case Excel workbook');
    if isequal(fileName,0), error('No workbook selected.'); end
    filePath = fullfile(folder,fileName);
end

T = readtable(filePath,'Sheet','All_40_cases','VariableNamingRule','preserve');
caseID = string(T.('Case'));
group = string(T.('Source test group'));
status = string(T.('Fit status'));
eD = T.('Fitted effective e/D');
Tco2 = T.('Spencer mean CO₂ temperature (°C)');
htc = T.('Spencer target HTC (W/m² K)');
hSmooth = T.('Numerical HTC at e/D=0 (W/m² K)');

% Excel blank cells may import as NaN; never treat unfitted cases as zero.
valid = contains(status,'HTC fitted within') & isfinite(eD) & eD>=0 & eD<=0.05;
groups = ["m12/q2","m12/q3","m15/q2","m15/q3"];
labels = {'m12 / q2','m12 / q3','m15 / q2','m15 / q3'};
colors = lines(4);

fprintf('Loaded %d total cases; %d valid fitted cases.\n',height(T),nnz(valid));

%% Figure 2: Four panels, each labeled by experimental case
figure('Color','w','Name','Individual Operating Groups');
for k=1:4
    idx = find(valid & group==groups(k));
    subplot(2,2,k); hold on;
    x = 1:numel(idx);
    scatter(x,eD(idx),65,colors(k,:),'filled');
    yline(mean(eD(idx)),'--','Mean','Color',colors(k,:),'LineWidth',1.3);
    xticks(x); xticklabels(cellstr(caseID(idx))); xtickangle(45);
    xlabel('Experimental case'); ylabel('Fitted e/D');
    title(sprintf('%s (n=%d)',labels{k},numel(idx)));
    ylim([0 0.055]); grid on; box on;
end

%% Figure 3: Group mean with standard deviation
means = nan(4,1);
sds = nan(4,1);
counts = zeros(4,1);

for k = 1:4
    v = eD(valid & group == groups(k));
    counts(k) = numel(v);
    means(k) = mean(v);
    sds(k) = std(v);
end

figure('Color','w','Name','Group Mean Roughness');
hold on;

for k = 1:4
    % Standard deviation error bars
    errorbar(k, means(k), sds(k), ...
        'Color', colors(k,:), ...
        'LineStyle', 'none', ...
        'LineWidth', 1.5, ...
        'CapSize', 12);

    % Mean roughness point
    scatter(k, means(k), 100, colors(k,:), ...
        'filled', 'MarkerEdgeColor', 'k');
end

xticks(1:4);
xticklabels(labels);
xlim([0.5 4.5]);

ylabel('Mean Fitted e/D');
title('Mean Fitted Roughness by Operating Condition (\pm1 SD)');

grid on;
box on;
hold off;

%% Figure 4: Roughness vs mean experimental CO2 temperature
figure('Color','w','Name','Roughness vs CO2 Temperature'); hold on;
for k=1:4
    idx = valid & group==groups(k) & isfinite(Tco2);
    scatter(Tco2(idx),eD(idx),70,colors(k,:),'filled', ...
        'DisplayName',labels{k});
end
xlabel('Spencer mean CO_2 temperature (°C)');
ylabel('Fitted effective relative roughness, e/D');
title('Fitted Roughness vs Mean CO_2 Temperature');
legend('Location','best'); grid on; box on; ylim([0 0.055]);

%% Figure 5: Fitted Roughness vs Experimental Target HTC
figure('Color','w','Name','Roughness vs Target HTC');
hold on;

for k = 1:4
    idx = valid & group == groups(k) & isfinite(htc);

    scatter(htc(idx),eD(idx),70,colors(k,:),'filled', ...
        'DisplayName',labels{k});
end

xlabel('Spencer Target HTC (W/m^2 K)');
ylabel('Fitted Effective Relative Roughness, e/D');
title('Fitted Roughness vs Experimental Target HTC');

legend('Location','best');
grid on;
box on;
ylim([0 0.055]);
hold off;

%% Figure 6: Old vs New Fitted Relative Roughness

% Locate old DH Excel file in Downloads folder
oldFilePath = fullfile(getenv('USERPROFILE'),'Downloads', ...
    'Spencer-vs-Numerical-HTC-All-40-Cases.xlsx');

% macOS/Linux fallback
if ~isfile(oldFilePath)
    oldFilePath = fullfile(getenv('HOME'),'Downloads', ...
        'Spencer-vs-Numerical-HTC-All-40-Cases.xlsx');
end

% File selection fallback
if ~isfile(oldFilePath)
    [fileName,folder] = uigetfile('*.xlsx', ...
        'Select the old DH Excel workbook');
    if isequal(fileName,0)
        error('No old DH workbook selected.');
    end
    oldFilePath = fullfile(folder,fileName);
end


% Read old DH data
Told = readtable(oldFilePath, ...
    'Sheet','All 40 cases', ...
    'VariableNamingRule','preserve');

oldCaseID = string(Told.('Case'));
oldED = Told.('Fitted effective e/D');

% Match old and new data using case names
[matched,oldIdx] = ismember(caseID,oldCaseID);

% Store matched old roughness values
oldEDmatched = nan(size(eD));
oldEDmatched(matched) = oldED(oldIdx(matched));

% Only plot cases with valid new e/D and available old e/D
compareValid = valid & matched & isfinite(oldEDmatched);

figure('Color','w','Name','Old vs New Relative Roughness');
hold on;

% Plot each operating condition
for k = 1:4
    idx = compareValid & group == groups(k);

    scatter(eD(idx),oldEDmatched(idx),75,colors(k,:), ...
        'filled', ...
        'DisplayName',labels{k});
end

% Equal scaling on both axes
axisMax = max([eD(compareValid); oldEDmatched(compareValid)]);
axisMax = ceil(axisMax*10)/10;

xlim([0 axisMax]);
ylim([0 axisMax]);
axis square;

% Perfect agreement reference line
plot([0 axisMax],[0 axisMax],'k--', ...
    'LineWidth',1.3, ...
    'DisplayName','Old e/D = New e/D');

xlabel('New Fitted Relative Roughness, e/D');
ylabel('Old DH Fitted Relative Roughness, e/D');
title('Old DH vs New Fitted Relative Roughness');

legend('Location','best');
grid on;
box on;
hold off;

%% Figures 7 and 8: Bulk-temperature plots matching the supplied PDF
% Figure 7: all 40 original fitted e/D values (red).
% Figure 8: revised fits only (blue), y-axis 0..0.01; off-scale counts labeled.
% Do NOT filter original fits using 'valid': that would omit original cases
% whose revised e/D could not be fitted.
%
% Tpc is the Cp maximum at fixed pressure, not the fixed critical temperature.
% HEOS references at the arithmetic-mean inlet pressure of ALL cases per class.
% Source: results/dimpled_trends/Pseudocritical_temperature_references.csv.
referencePressureMPa = [8.048454545454545, 8.0639, ...
                        8.015555555555556, 8.0473];
pseudoCriticalC = [34.95584910688052, 35.04559676121892, ...
                   34.764202319169044, 34.949131746483374];

pressureMPa = T.('Experimental inlet pressure (MPa)');
reportedMassFlux = T.('Reported CO₂ mass flux (kg/m² s)');
nominalFlux = T.('Nominal experimental heat-flux group (W/m²)');

% Refuse to silently apply the saved Tpc references to a different dataset.
for k = 1:4
    rows = group == groups(k) & isfinite(pressureMPa);
    if ~any(rows)
        error('No finite inlet pressures available for group %s.',char(groups(k)));
    end
    actualMeanMPa = mean(pressureMPa(rows));
    if abs(actualMeanMPa-referencePressureMPa(k)) > 5e-5
        error(['Group %s has mean inlet pressure %.6f MPa, but its saved ' ...
            'Tpc reference uses %.6f MPa. Recalculate Tpc for this dataset.'], ...
            char(groups(k)),actualMeanMPa,referencePressureMPa(k));
    end
end

oldValid = matched & isfinite(oldEDmatched) & isfinite(Tco2);
newTempValid = valid & isfinite(Tco2);
plotOutputDir = fullfile(fileparts(filePath),'Dimpled_Temperature_Plots');
if ~isfolder(plotOutputDir), mkdir(plotOutputDir); end

for overlayMode = 0:1
    if overlayMode == 0
        figureName = 'Original 40 Fits vs Bulk Temperature';
        subtitleText = 'Original 40 DH-based fitted values';
        outputStem = 'eD_vs_Bulk_Temperature_Original_40';
    else
        figureName = 'Corrected Fits Only vs Bulk Temperature';
        subtitleText = ['Revised smooth-corrected fits only (blue; 24 matched cases); ' ...
            'y-axis zoom: 0 to 0.01'];
        outputStem = 'eD_vs_Bulk_Temperature_Corrected_Blue_0_to_001';
    end
    fig = figure('Color','w','Name',figureName, ...
        'Position',[100 100 1400 850]);

    for k = 1:4
        ax = subplot(2,2,k); hold(ax,'on');
        idxOld = oldValid & group == groups(k);
        idxNew = newTempValid & group == groups(k);
        if overlayMode == 0
            scatter(ax,Tco2(idxOld),oldEDmatched(idxOld),65,'r','filled', ...
                'DisplayName',sprintf('Original fits (n=%d)',nnz(idxOld)));
        end
        if overlayMode == 1
            scatter(ax,Tco2(idxNew),eD(idxNew),60,'b','filled', ...
                'DisplayName',sprintf('Corrected fits (n=%d)',nnz(idxNew)));
        end

        if overlayMode == 0
            yline(ax,0.05,'--','Color',[0.4 0.4 0.4],'LineWidth',1.1, ...
                'DisplayName','Reference e/D = 0.05');
        else
            above = nnz(idxNew & eD > 0.01);
            text(ax,0.98,0.05,sprintf('%d of %d blue values above 0.01 (off-scale)', ...
                above,nnz(idxNew)),'Units','normalized', ...
                'HorizontalAlignment','right','FontSize',8);
        end
        xline(ax,pseudoCriticalC(k),'--','Color',[0 0.51 0.29], ...
            'LineWidth',1.8, ...
            'DisplayName',sprintf('T_{pc} = %.2f C at %.3f MPa', ...
                pseudoCriticalC(k),referencePressureMPa(k)));

        idxGroup = group == groups(k);
        finiteG = reportedMassFlux(idxGroup & isfinite(reportedMassFlux));
        finiteQ = nominalFlux(idxGroup & isfinite(nominalFlux));
        title(ax,sprintf('%s | G = %.1f-%.1f kg/(m^2 s) | nominal flux = %.0f kW/m^2', ...
            char(groups(k)),min(finiteG),max(finiteG),mean(finiteQ)/1000), ...
            'FontSize',10);
        xlabel(ax,'Reported mean CO_2 bulk temperature (°C)');
        ylabel(ax,'Fitted effective e/D');
        xlim(ax,[25 39.5]); xticks(ax,26:2:38);
        if overlayMode == 0
            ylim(ax,[0 0.70]); yticks(ax,0:0.1:0.7);
        else
            ylim(ax,[0 0.01]); yticks(ax,0:0.002:0.01);
        end
        grid(ax,'on'); box(ax,'on');
        legend(ax,'Location','northwest','FontSize',8);
        hold(ax,'off');
    end

    sgtitle({'Fitted e/D vs bulk CO_2 temperature by operating-condition class', ...
        subtitleText},'FontSize',14);
    noteText = ['Green line: CO_2 Cp maximum at each class mean inlet pressure (HEOS). ' ...
        'Individual runs have pressure-dependent Tpc. e/D is fitted, not measured roughness.'];
    if overlayMode == 1
        noteText = [noteText ' Values above 0.01 are off-scale; 16 revised cases are unfitted.'];
    end
    annotation(fig,'textbox',[0.07 0.005 0.86 0.045], ...
        'String',noteText,'EdgeColor','none','HorizontalAlignment','center', ...
        'FontSize',8,'Interpreter','none');

    exportgraphics(fig,fullfile(plotOutputDir,[outputStem '.pdf']), ...
        'ContentType','vector');
    exportgraphics(fig,fullfile(plotOutputDir,[outputStem '.png']), ...
        'Resolution',300);
end

fprintf('\nSaved the two temperature plots as PDF and PNG in:\n%s\n',plotOutputDir);
