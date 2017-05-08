function plotConverTrend( FEsEachGen, bestFit, saveFigPath, figName )
% Plot the convergence trend.
%   Parameters:
%   FEsEachGen          - The number of fitness evaluations of each generation
%                       [row vector]
%   bestFit             - The best fitness values, can be bestFitSoFar or bestFitEachGen
%                       [row vector]
%   saveFigPath         - The path to save figure
%                       [positive scalar]
%   figName             - The name of the saved figure
%                       [string]


figure('Visible', 'off');
plot(FEsEachGen, bestFit, 'b');
xlabel('FEs');
ylabel('f(x)');
grid on;
print([saveFigPath, filesep, figName], '-depsc');
close;

end

