function GA_core( run_num, func_names, func_num, options, optionalArgs )
% Perform the core part of Genetic Algorithm (GA).
%   Parameters:
%   run_num             - Run number
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


% set different rand number generator in each run
rng(round(sum(clock) * 1000), 'twister');

if isfield(optionalArgs, 'fst_num')
    % Load variables of the first generation, struct array
    variables = optionalArgs.variables;
    FEs = variables.FEs;
    gens = variables.gens;
    if isfield(variables, 'params')
        params = variables.params; % exogenous params struct
        % params is also a struct array
        if ~isstruct(params)
            error('params should be a struct array!');
        end
    end
else
    % Init several local variables
    FEs = 0; % total fitness evaluations
    gens = 0; % total number of generations
    bestFitEachGen = zeros(1); % init best fitness values of each generation
    bestXEachGen = cell(1); % init best x of each generation
    FEsEachGen = zeros(1); % init number of fitness evaluations of each generation
    bestXSoFar = cell(1); % init best x so far
    bestFitSoFar = zeros(1); % init best fitness values so far
end


% Init the first generation
[x, fit, saveGenPath] = initFstGen(optionalArgs, options, func_names, func_num, run_num);

if ~isfield(optionalArgs, 'fst_num')
    % find the best fitness value and its index
    [minFit, bestIdx] = min(fit);
    bestFit = minFit;
    % find the best individual
    bestX = x(bestIdx, :);
    % update local variables
    FEs = FEs + options.PopulationSize;
    gens = gens + 1;
    bestFitEachGen(gens) = bestFit;
    bestFitSoFar(gens) = bestFit;
    bestXEachGen{gens} = bestX;
    bestXSoFar{gens} = bestX;
    FEsEachGen(gens) = FEs;
    % build a struct array to save inner variables
    variables.FEs = FEs;
    variables.gens = gens;
    % save this generation to file
    saveEachGen(gens, saveGenPath, x, fit, variables);
end

% load exogenous adaptive parameters
% no exogenous parameters for GA

% Perform the optimization process
if isfield(optionalArgs, 'fst_num')
    % just generate one offspring if provided the field 'fst_num'
    FEs = options.MaxFEs - 1;
end

% do the main loop
while FEs < options.MaxFEs
    % calculate how many different kinds of offspring will be generated
    nEliteKids = options.EliteCount; % number of elite kids
    % number of crossover kids
    nXoverKids = round(options.CrossoverFraction * (options.PopulationSize - nEliteKids));
    % calculate how many parents need be selected to perform crossover
    % crossover: two parents generate two kids
    % need to make sure nParentsXover is even
    if ~rem(nXoverKids, 2)
        nParentsXover = nXoverKids;
    else
        nParentsXover = nXoverKids - 1;
    end
    nXoverKids = nParentsXover; % update the new number of crossover kids
    % number of mutation kids
    nMutationKids = options.PopulationSize - nEliteKids - nXoverKids;
    
    % select parents for crossover
    % fitness scaling, default is 'rank'
    fitScaled = fitScaling(fit, options.PopulationSize, options);
    % parent selection, default is 'SUS'
    % the returned parents are indices of selected parents
    parents = parentSelection(fitScaled, nXoverKids, options);
    % randomly shuffle the selected parents' indices
    parents = parents(randperm(length(parents)));
    
    % the crossover kids
    % perform crossover: two parents generate two kids
    xoverKids = crossover(x, parents, options);
    
    % mutation on xoverKids: one parent generate one kid
    mutationIndices = randperm(nXoverKids, nMutationKids);
    mutationKids = mutation(xoverKids(mutationIndices, :), options);
    
    % the elite kids
    [~, ind] = sort(fit);
    elitsKids = x(ind(1:nEliteKids), :);
    fitElite = fit(ind(1:nEliteKids));
    
    % survival selection using the generational replacement to form the new x
    xEval = [xoverKids; mutationKids];
    
    % bound constraint (default: absorb)
    xEval = boundCons(xEval, options);
    
    % evaluate the new population
    fitEval = benchmark_func(xEval, func_num);
    
    % add the elite kids
    x = [elitsKids; xEval];
    fit = [fitElite; fitEval];
    
    % find the best fitness value and its index
    [minFit, bestIdx] = min(fit);
    bestFit = minFit;
    % find the best individual
    bestX = x(bestIdx, :);
    
    if isfield(optionalArgs, 'fst_num')
        % get the best fitness value of the neighbour
        f_ij_b = bestFit;
        
        FEs = FEs + options.PopulationSize - nEliteKids;
        gens = gens + 1;
    else
        % update local variables
        FEs = FEs + options.PopulationSize - nEliteKids;
        gens = gens + 1;
        bestFitEachGen(gens) = bestFit;
        bestXEachGen{gens} = bestX;
        % judge best fitness value so far
        if bestFit < bestFitSoFar(gens-1)
            bestFitSoFar(gens) = bestFit;
            bestXSoFar{gens} = bestX;
        else
            bestFitSoFar(gens) = bestFitSoFar(gens-1);
            bestXSoFar(gens) = bestXSoFar(gens-1);
        end
        FEsEachGen(gens) = FEs;
        % build a struct array to save inner variables
        variables.FEs = FEs;
        variables.gens = gens;
        % save this generation to file
        saveEachGen(gens, saveGenPath, x, fit, variables);
        
        % ignore the last generation if: leftFEs < populationSize
        if options.MaxFEs - FEs < options.PopulationSize - nEliteKids
            break;
        end
        
    end
    
end

if isfield(optionalArgs, 'fst_num')
    % save f_ij_b to a .mat file
    if ~isnumeric(f_ij_b)
        fprintf('The f_ij_b should be a numeric value!');
    end
    save([optionalArgs.saveImpPath, filesep, num2str(optionalArgs.fst_num), ...
        '_', num2str(run_num), '.mat'], 'f_ij_b');
else
    % save variables to .mat file and plot the convergence curve
    bestFitness(bestFitEachGen, bestFitSoFar, bestXEachGen, bestXSoFar, FEsEachGen, ...
        func_names, func_num, run_num, options);
end

end

