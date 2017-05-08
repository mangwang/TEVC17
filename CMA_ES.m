function CMA_ES( run_num, func_names, func_num, options, optionalArgs )
% Framework for CMA Evolution Strategy (CMA-ES).
%   Parameters:
%   run_num            - Run number
%                       [positive scalar]
%   func_names          - Function names
%                       [cell array of strings]
%   func_num            - The number of optimization functions
%                       [positive scalar]
%   options             - The options set by setOptions()
%                       [struct array]
%   optionalArgs        - The optional arguments, feasible fields include
%                         'InitPopulation', 'InitFitness', 'fst_num', 'variables'
%                       [struct array]
%   Fields of optionalArgs:
%   InitPopulation      - The initial population used in seeding the GA
%                         algorithm
%                       [ Matrix | [] ]
%   InitFitness         - The initial fitness values used to determine fitness
%                       [ column vector | [] ]
%   fst_num             - The number of first generation, it is also a flag
%                         to check whether to generate just one offspring
%                       [positive scalar]
%   variables           - Inner variables provided by the first generation
%                       [struct array]


if isfield(optionalArgs, 'fst_num') && run_num == 1
    % load the path of global best fitness value in the history
    gBestPath = ['result', filesep, 'conver_trend', filesep, 'dim_', num2str(options.Dim),...
        filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
        'run_', num2str(optionalArgs.run_num_normal)];
    % read the vector of bestFitSoFar
    gBestStr = load([gBestPath, filesep, 'bestFitSoFar.mat']);
    gBest = gBestStr.bestFitSoFar;
    optionalArgs.gBest = gBest; % save to optionalArgs
    
    % get the best fitness value of the neighbour
    f_ij_b = gBest(optionalArgs.fst_num + 1);
    
    % save f_ij_b to a .mat file
    if ~isnumeric(f_ij_b)
        fprintf('The f_ij_b should be a numeric value!');
    end
    save([optionalArgs.saveImpPath, filesep, num2str(optionalArgs.fst_num), '_1.mat'], ...
        'f_ij_b');
    
end

% run the core
if ~isfield(optionalArgs, 'nNeighbours') || (isfield(optionalArgs, 'nNeighbours') && run_num > 1)
    CMA_ES_core(run_num, func_names, func_num, options, optionalArgs);
end

end

