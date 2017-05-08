function [ x, fit, saveGenPath ] = initFstGen( optionalArgs, options, func_names, ...
    func_num, run_num )
% Initialise the first generation including genotype and phenotype..
%   Parameters:
%   optionalArgs        - The optional arguments, feasible fields include
%                         all fields in options and 'InitPopulation',
%                         'InitFitness', 'fst_num'
%                       [struct array | See definitions of EA_setOptions()]
%   options             - The options get from EA_setOptions() function
%                       [struct array | See definitions of EA_setOptions()]
%   func_name           - Function names
%                       [cell array of strings]
%   func_num            - The number of optimization functions
%                       [positive scalar]
%   run_num             - The number of run times
%                       [positive scalar]
%   Fields of optionalArgs:
%   InitPopulation      - The initial population used in seeding the GA
%                         algorithm
%                       [ Matrix | [] ]
%   InitFitness         - The initial fitness values used to determine fitness
%                       [ column vector | [] ]
%   fst_num             - The number of first generation
%                       [positive scalar]


if isfield(optionalArgs, 'InitPopulation') && isfield(optionalArgs, 'InitFitness')
    % load x and fitness
    x = optionalArgs.InitPopulation;
    fit = optionalArgs.InitFitness;
else
    % generate first generation that randomly and uniformly distributed in the feasible region
    lb = options.PopInitRange(1);
    ub = options.PopInitRange(2);
    x = lb + (ub - lb) * rand(options.PopulationSize, options.Dim);
    fit = benchmark_func(x, func_num);
end

% save path
saveGenPath = ['result', filesep, 'raw_data', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
    'run_', num2str(run_num)];

end

