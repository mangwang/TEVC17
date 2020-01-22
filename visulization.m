function visulization( run_num_normal, func_names, func_num, options )
% Show several figures for the relations of indicators and performance.
%   Parameters:
%   runNumNormal        - The number of normal runs
%                       [positive scalar]
%   func_names          - Function names
%                       [cell array of strings]
%   func_num            - The number of optimization functions
%                       [positive scalar]
%   options             - The options set by setOptions()
%                       [struct array]


% load all the indicators
loadEvoPath = ['result', filesep, 'evolvability', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
    'run_', num2str(run_num_normal)];
indStr = load([loadEvoPath, filesep, 'indicators.mat']);
eppAllGens = indStr.eppAllGens;  % epp
eapAllGens = indStr.eapAllGens;  % eap
evpAllGens = indStr.evpAllGens;  % evp

nCalGens = length(evpAllGens);  % obtain the number of generations for calculating measures

% path to save the figure
saveVisPath = ['result', filesep, 'visulization', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}];
if ~isdir(saveVisPath)
    mkdir(saveVisPath);
end

% get the benchmark results
readConvBase = ['result', filesep, 'benchmark', filesep, 'conver_trend_all_runs', filesep, ...
    'dim_', num2str(options.Dim), filesep, func_names{func_num}, filesep, options.AlgoName];
% read the variables from bestFitAllRuns
convStr = load([readConvBase, filesep, 'bestFitAllRuns.mat']);
bestFitMed = convStr.bestFitMed;  % median
bestFitAvg = convStr.bestFitAvg;  % mean
FEsEachGen = convStr.FEsEachGen;  % FEs

% set global optimal to 0
if func_num == 14
    bestFitAvgPlot = bestFitAvg(1:(nCalGens+1));
    bestFitMedPlot = bestFitMed(1:(nCalGens+1));
else
    bestFitAvgPlot = bestFitAvg(1:(nCalGens+1)) + 600 - 100 * (func_num - 1);
    bestFitMedPlot = bestFitMed(1:(nCalGens+1)) + 600 - 100 * (func_num - 1);
end
avgIdx = find(bestFitAvgPlot == bestFitAvgPlot(end), 1);
medIdx = find(bestFitMedPlot == bestFitMedPlot(end), 1);

nGens = min([nCalGens, avgIdx, medIdx]);  % the number of generations for ploting

% plot four subfigures: epp + eap + evp + conv_trend
saveAlgoName = strrep(options.AlgoName, '_', '-');
fig = figure('Visible', 'off');

ax1 = subplot(4, 1, 1);
semilogy(ax1, FEsEachGen(1:nGens), eppAllGens(1:nGens));
set(ax1, 'XTickLabel', {' '}, 'TickLength', [0.001, 0.001]);
ylim(ax1, [0, 1.1]);
ylabel('$log(epp)$', 'interpreter', 'latex');
title(saveAlgoName);
grid on;

ax2 = subplot(4, 1, 2);
semilogy(ax2, FEsEachGen(1:nGens), eapAllGens(1:nGens));
set(ax2, 'XTickLabel', {' '}, 'TickLength', [0.001, 0.001]);
ylabel('$log(eap)$', 'interpreter', 'latex');
grid on;

ax3 = subplot(4, 1, 3);
semilogy(ax3, FEsEachGen(1:nGens), evpAllGens(1:nGens));
set(ax3, 'XTickLabel', {' '}, 'TickLength', [0.001, 0.001]);
ylabel('$log(evp)$', 'interpreter', 'latex');
grid on;

ax4 = subplot(4, 1, 4);
% plot mean
semilogy(ax4, FEsEachGen(1:nGens), bestFitAvgPlot(1:nGens));
hold on;
% plot median
semilogy(ax4, FEsEachGen(1:nGens), bestFitMedPlot(1:nGens));
legend({'mean', 'median'}, 'FontSize', 6, 'Location', 'best');
xlabel('The Number of FEs');
ylabel('$log(f-f^*)$', 'interpreter', 'latex');
grid on;

% save figure
if strcmp(options.form, 'eps')
    print(fig, [saveVisPath, filesep, saveAlgoName], '-depsc', '-r1200');
else
    print(fig, [saveVisPath, filesep, saveAlgoName], '-dpdf', '-r1200');
end
close(fig);

end

