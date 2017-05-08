function CoDE_core( run_num, func_names, func_num, options, optionalArgs )
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


% ------------------------------------------------------------------------------------------
% Please refer to the official website:
% http://dces.essex.ac.uk/staff/zhang/code/codealgorithm/codealgorithm.rar
% ------------------------------------------------------------------------------------------


% set different rand number generator in each run
rng(round(sum(clock) * 1000), 'twister');

% fixed parameters
NP = options.PopulationSize;
D = options.Dim;
nStras = 3;  % the number of strategies
F = [1, 1, 0.8];
CR = [0.1, 0.9, 0.2];

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
% no exogenous parameters for CoDE

% Perform the optimization process
if isfield(optionalArgs, 'fst_num')
    % just generate one offspring if provided the field 'fst_num'
    FEs = options.MaxFEs - 1;
end

% do the main loop
while (FEs < options.MaxFEs)
    % generate mutaMatrix for mutation
    newX = zeros(NP * nStras, D);
    newFit = zeros(NP, nStras);
    randp = randperm(nStras);
    for i = 1:nStras
        selFCR = randp(i);
        % generate crossoverIndex
        selIndex = rand(NP, D) < CR(selFCR);
        randI = randi(D, NP, 1);
        for k = 1:NP
            selIndex(k, randI(k)) = true;
        end
        sel = ~selIndex;
        if i == 1
            % rand/1/bin
            nMuta = 3;
            mutaMatrix = zeros(nMuta, NP);
            for j = 1:NP
                mutaMatrix(:, j) = randperm(NP, nMuta);
            end
            tempX = x(mutaMatrix(1, :), :) + F(selFCR) .* (x(mutaMatrix(2, :), :) - x(mutaMatrix(3, :), :));
            tempX(sel) = x(sel);
        elseif i == 2
            % rand/2/bin
            nMuta = 5;
            mutaMatrix = zeros(nMuta, NP);
            for j = 1:NP
                mutaMatrix(:, j) = randperm(NP, nMuta);
            end
            tempX = x(mutaMatrix(1, :), :) + repmat(rand(NP, 1), 1, D) .* (x(mutaMatrix(2, :), :) - x(mutaMatrix(3, :), :))...
                + F(selFCR) .* (x(mutaMatrix(4, :), :) - x(mutaMatrix(5, :), :));
            tempX(sel) = x(sel);
        elseif i == 3
            % current to rand/1
            mutaMatrix = randi(NP, 3, NP);
            tempX = x + repmat(rand(NP, 1), 1, D) .* (x(mutaMatrix(1, :), :) - x) + ...
                F(selFCR) .* (x(mutaMatrix(2, :), :) - x(mutaMatrix(3, :), :));
        end
        newX((1:NP) + NP * (i - 1), :) = tempX;
    end
    
    % bound constraint  (default: absorb)
    newX = boundCons(newX, options);
    
    newFit(:) = benchmark_func(newX, func_num);
    
    % select the best strategies
    [selFit, minInd] = min(newFit, [], 2);
    minInd = (minInd - 1) * NP + (1:NP)';
    selX = newX(minInd, :);
    minInd = selFit <= fit;
    % update x and fit
    x(minInd, :) = selX(minInd, :);
    fit(minInd) = selFit(minInd);
    
    % find the best fitness value and its index
    [minFit, bestIdx] = min(fit);
    bestFit = minFit;
    % find the best individual
    bestX = x(bestIdx, :);
    
    if isfield(optionalArgs, 'fst_num')
        % get the best fitness value of the neighbour
        f_ij_b = bestFit;
        
        FEs = FEs + 3 * NP;  % 3*NP
        gens = gens + 1;
    else
        % update local variables
        FEs = FEs + 3 * NP;
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
        if options.MaxFEs - FEs < 3 * options.PopulationSize
            break;
        end
        
    end
    
end

if isfield(optionalArgs, 'fst_num')
    if ~isnumeric(f_ij_b)
        fprintf('The f_ij_b should be a numeric value!');
    end
    % save f_ij_b to a .mat file
    save([optionalArgs.saveImpPath, filesep, num2str(optionalArgs.fst_num), ...
        '_', num2str(run_num), '.mat'], 'f_ij_b');
else
    % save variables to .mat file and plot the convergence curve
    bestFitness(bestFitEachGen, bestFitSoFar, bestXEachGen, bestXSoFar, FEsEachGen, ...
        func_names, func_num, run_num, options);
end

end
