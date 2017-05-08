function ABC_core( run_num, func_names, func_num, options, optionalArgs )
% Perform the core part of Artificial Bee Colony (ABC).
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
% http://mf.erciyes.edu.tr/abc/
% ------------------------------------------------------------------------------------------


% set different rand number generator in each run
rng(round(sum(clock) * 1000), 'twister');

% fixed parameters
D = options.Dim;
lb = options.PopInitRange(1);
ub = options.PopInitRange(2);
foodNumber = options.PopulationSize;  % The number of food sources equals the half of the colony size
limit = 100;    % A food source which could not be improved through "limit" trials is abandoned by its employed bee

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
    % find the best food source
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
% no exogenous parameters for ABC

% Perform the optimization process
if isfield(optionalArgs, 'fst_num')
    % just generate one offspring if provided the field 'fst_num'
    FEs = options.MaxFEs - 1;
end

fitness = calculateFitness(fit);  % fitness of ABC
trial = zeros(1, foodNumber);  % reset trial counters
while FEs < options.MaxFEs
    % ---- EMPLOYED BEE PHASE ---- %
    for i = 1:foodNumber    % the number of employed bee is same as the food number
        % The parameter to be changed is determined randomly
        param2Change = fix(rand*D)+1;        % fix function
        % A randomly chosen solution is used in producing a mutant solution of the solution i
        neighbour = fix(rand*(foodNumber))+1;
        % Randomly selected solution must be different from the solution i
        while neighbour == i
            neighbour = fix(rand*(foodNumber))+1;
        end
        sol = x(i,:);
        %  v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})
        sol(param2Change) = x(i,param2Change)+(x(i,param2Change)-x(neighbour,param2Change))*(rand-0.5)*2;
        
        %  if generated parameter value is out of boundaries, it is shifted onto the boundaries
        sol = boundCons(sol, options);
        
        % evaluate new solution
        objValSol = benchmark_func(sol, func_num);
        FEs = FEs + 1;
        fitnessSol = calculateFitness(objValSol);
        
        % a greedy selection is applied between the current solution i and its mutant
        if fitnessSol > fitness(i) % If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i
            x(i,:) = sol;
            fitness(i) = fitnessSol;
            fit(i) = objValSol;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1; % if the solution i can not be improved, increase its trial counter
        end
    end
    % ---- CalculateProbabilities ---- %
    % A food source is chosen with the probability which is proportioal to its quality
    % Different schemes can be used to calculate the probability values
    % For example prob(i)=fitness(i)/sum(fitness)
    % or in a way used in the metot below prob(i) = a*fitness(i)/max(fitness)+b
    % probability values are calculated by using fitness values and normalized by dividing maximum fitness value
    prob = 0.9 .* fitness ./ max(fitness) + 0.1;  % the probability of onlooker bees to select foods
    % ---- ONLOOKER BEE PHASE ---- %
    i = 1;
    t = 0;
    while t < foodNumber
        if rand < prob(i)
            t = t + 1;
            % The parameter to be changed is determined randomly
            param2Change = fix(rand*D)+1;
            % A randomly chosen solution is used in producing a mutant solution of the solution i
            neighbour = fix(rand*(foodNumber))+1;
            % Randomly selected solution must be different from the solution i
            while(neighbour == i)
                neighbour = fix(rand*(foodNumber))+1;
            end
            sol = x(i,:);
            % v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})
            sol(param2Change) = x(i,param2Change)+(x(i,param2Change)-x(neighbour,param2Change))*(rand-0.5)*2;
            % if generated parameter value is out of boundaries, it is shifted onto the boundaries
            sol = boundCons(sol, options);
            % evaluate new solution
            objValSol = benchmark_func(sol, func_num);
            FEs = FEs + 1;
            fitnessSol = calculateFitness(objValSol);
            % a greedy selection is applied between the current solution i and its mutant
            if fitnessSol > fitness(i) % If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i
                x(i,:) = sol;
                fitness(i) = fitnessSol;
                fit(i) = objValSol;
                trial(i) = 0;
            else
                trial(i) = trial(i)+1; % if the solution i can not be improved, increase its trial counter
            end
        end
        i = i+1;
        if i == foodNumber+1
            i = 1;
        end
    end
    % ---- SCOUT BEE PHASE ---- %
    % determine the food sources whose trial counter exceeds the "limit" value.
    % In Basic ABC, only one scout is allowed to occur in each cycle
    ind = find(trial==max(trial));
    ind = ind(end);
    if trial(ind) > limit
        trial(ind) = 0;
        sol = rand(1, D) .* (ub - lb) + lb;
        objValSol = benchmark_func(sol, func_num);
%         FEs = FEs+1;
        fitnessSol = calculateFitness(objValSol);
        x(ind,:) = sol;
        fitness(ind) = fitnessSol;
        fit(ind) = objValSol;
    end
    
    % find the best fitness value and its index
    [minFit, bestIdx] = min(fit);
    bestFit = minFit;
    % find the best individual
    bestX = x(bestIdx, :);
    
    if isfield(optionalArgs, 'fst_num')
        % get the best fitness value of the neighbour
        f_ij_b = bestFit;
        
        gens = gens + 1;
    else
        % update local variables
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
        if options.MaxFEs - FEs < 2 * foodNumber
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


function ABC_Fit = calculateFitness(fit)
% calculate fitness value in ABC algorithm.
ABC_Fit = zeros(size(fit));
ind = find(fit >= 0);
ABC_Fit(ind) = 1 ./ (fit(ind) + 1);
ind = find(fit < 0);
ABC_Fit(ind) = 1 + abs(fit(ind));
end
