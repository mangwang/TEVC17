function SPSO2011_core(run_num, func_names, func_num, options, optionalArgs)
% Perform the core part of Standard PSO 2011 (SPSO2011).
%   Parameters:
%   run_num             - The number of run times
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
% http://www.particleswarm.info/SPSO2011_matlab.zip
% ------------------------------------------------------------------------------------------


% set different rand number generator in each run
rng(round(sum(clock) * 1000), 'twister');

% fixed parameters
lb = options.PopInitRange(1);  % lower bound
ub = options.PopInitRange(2);  % upper bound
NP = options.PopulationSize;  % population size
D = options.Dim;  % dimension of optimization function
w = 1 / (2 * log(2));
c1 = 0.5 + log(2);
c2 = c1;
k = 3;
p = 1 - (1 - 1 / NP) ^ k;  % probability to be an informant

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
    FEs = FEs + NP;
    gens = gens + 1;
    bestFitEachGen(gens) = bestFit;
    bestFitSoFar(gens) = bestFit;
    bestXEachGen{gens} = bestX;
    bestXSoFar{gens} = bestX;
    FEsEachGen(gens) = FEs;
    % build a struct array to save inner variables
    variables.FEs = FEs;
    variables.gens = gens;
    % init particles' velocity as a uniformly distributed random number in the
    % range ((lb-x), (ub-x))
    vel = (lb - x) + (ub - lb) * rand(NP, D);
    % init personal best (pbest) and fitness values
    pbest = x;
    pFit = fit;
    % init global best (gbest) and fitness values
    gbest = bestX;
    gFit = minFit;
    % init nighs for the randomize topology
    nighs = rand(NP) < p;
    for s = 1:NP
        nighs(s, s) = true;
    end
    % init parameters for the first generation
    params = struct('vel', vel, 'pbest', pbest, 'pFit', pFit,...
        'gbest', gbest, 'gFit', gFit, 'nighs', nighs);
    variables.params = params;
    % save this generation to file
    saveEachGen(gens, saveGenPath, x, fit, variables);
end


% load exogenous adaptive parameters
if isfield(variables, 'params')
    params = variables.params; % exogenous adaptive params struct
    % params is also a struct array
    if ~isstruct(params)
        error('params should be a struct array!');
    end
    % load parameters first
    vel = params.vel;
    pbest = params.pbest;
    pFit = params.pFit;
    gbest = params.gbest;
    gFit = params.gFit;
    nighs = params.nighs;
end

% Perform the optimization process
if isfield(optionalArgs, 'fst_num')
    % just generate one offspring if provided the field 'fst_num'
    FEs = options.MaxFEs - 1;
end

while (FEs < options.MaxFEs)
    % find the best informant
    nighsFit = repmat(pFit', NP, 1);
    nighsFit(~nighs) = NaN;
    [~, gbestInd] = min(nighsFit, [], 2);
    
    % define the gravity center G
    sw = 3 * ones(NP, 1);
    c22 = c2 * ones(NP, 1);
    for i = 1:NP
        if gbestInd(i) == i
            sw(i) = 2;
            c22(i) = 0;
        end
    end
    G = x + (c1 * (pbest - x) + repmat(c22, 1, D) .* (pbest(gbestInd, :) - x)) ./ repmat(sw, 1, D);
    
    % generate random points in the hypersphere aound G (uniform distribution)
    radius = sqrt(sum((G - x) .^ 2, 2));  % radius is the Euclidean norm of (G-x)
    randSphere = randn(NP, D);
    normRandSphere = sqrt(sum(randSphere .^ 2, 2));
    newPop = randSphere ./ repmat(normRandSphere, 1, D) .* repmat(rand(NP, 1) .* radius, 1, D) + G;
    
    % update the velocity = w * vel + random_points - x
    vel = w * vel + newPop - x;
    
    % update x
    x = x + vel;
    
    % bound constraint: absorb also update vel
    ind = x < lb;
    x(ind) = lb;
    vel(ind) = -0.5 * vel(ind);
    ind = x > ub;
    x(ind) = ub;
    vel(ind) = -0.5 * vel(ind);
    
    fit = benchmark_func(x, func_num);
    
    % update personal best
    ind = fit < pFit;
    pFit(ind, :) = fit(ind, :);
    pbest(ind, :) = x(ind, :);
    
    % update global best
    [minFit, bestIdx] = min(fit);
    if minFit < gFit
        gFit = minFit;
        gbest = x(bestIdx, :);
        adaptNighs = false;
    else
        adaptNighs = true;
    end
    bestFit = gFit;
    bestX = gbest;
    
    % randomize topolpgy if adaptNighs is true
    if adaptNighs == true  % no improvement in the best solution
        nighs = rand(NP) < p;
        for s = 1:NP
            nighs(s, s) = true;
        end
    end
    
    if isfield(optionalArgs, 'fst_num')
        % get the best fitness value of the neighbour
        f_ij_b = bestFit;
        
        FEs = FEs + NP;
        gens = gens + 1;
    else
        % update local variables
        FEs = FEs + NP;
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
        params = struct('vel', vel, 'pbest', pbest, 'pFit', pFit,...
            'gbest', gbest, 'gFit', gFit, 'nighs', nighs);
        variables.params = params; % save params in the first generation
        % save this generation to file
        saveEachGen(gens, saveGenPath, x, fit, variables);
        
        % ignore the last generation if: leftFEs < populationSize
        if options.MaxFEs - FEs < NP
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

