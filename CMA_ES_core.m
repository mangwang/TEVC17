function CMA_ES_core( run_num, func_names, func_num, options, optionalArgs )
% Perform the core part of CMA Evolution Strategy (CMA_ES).
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
% https://www.lri.fr/~hansen/purecmaes.m
% ------------------------------------------------------------------------------------------


% set different rand number generator in each run
rng(round(sum(clock) * 1000), 'twister');

% fixed parameters
lb = options.PopInitRange(1);  % lower bound
ub = options.PopInitRange(2);  % upper bound
lambda = options.PopulationSize;  % population size, offspring number
mu = lambda / 2;  % number of parents/points for recombination
mu = floor(mu);
n = options.Dim;  % dimension of optimization function

% Strategy parameter setting: Recombination
weights = log(mu+1/2) - log(1:mu);  % 1 * mu array for weighted recombination
weights = weights / sum(weights);  % normalize recombination weights array
mueff = sum(weights)^2 / sum(weights.^2);  % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4 + mueff/n) / (n + 4 + 2*mueff/n);  % time constant for cumulation for C
cs = (mueff + 2) / (n + mueff + 5);  % t-const for cumulation for sigma control
c1 = 2 / ((n + 1.3) ^ 2 + mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((n+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(n+1))-1) + cs;  % damping for sigma, usually close to 1

chiN = n^0.5*(1-1/(4*n)+1/(21*n^2));  % expectation of ||n(0,I)|| == norm(randn(n,1))

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
    % sort the fitness values
    [~, ind] = sort(fit);
    % find the index of best individual
    bestIdx = ind(1);
    bestFit = fit(bestIdx);  % the best fitness value
    bestX = x(bestIdx, :);  % the best individual
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
    % init exogenous adaptive parameters (#8)
    pc = zeros(1, n);  % evolution paths for C
    ps = zeros(1, n);   % evolution paths for sigma
    B = eye(n, n);  % B defines the coordinate system
    D = ones(n, 1);  % diagonal D defines the scaling
    C = B * diag(D.^2) * B';  % covariance matrix C
    invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
    xmean = weights * x(ind(1:mu), :);  % init xmean
    sigma = 0.3 * (ub - lb);  % coordinate wise standard deviation (step size)
    
    params = struct('B', B, 'D', D, 'xmean', xmean, 'sigma', sigma);
    variables.params = params; % save params in the first generation
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
    B = params.B;  % B defines the coordinate system
    D = params.D;  % diagonal D defines the scaling
    xmean = params.xmean;  % xmean
    sigma = params.sigma;  % coordinate wise standard deviation (step size)
end

% Perform the optimization process
if isfield(optionalArgs, 'fst_num')
    % just generate one offspring if provided the field 'fst_num'
    FEs = options.MaxFEs - 1;
end

while (FEs < options.MaxFEs)
    % mutation
    x = repmat(xmean, lambda, 1) + ...
        sigma * randn(lambda, n) * diag(D) * B';
    
    % bound constraint: resampling
    i = 1;
    while i <= lambda
        if any(x(i, :) < lb) || any(x(i, :) > ub)
            % resampling
            x(i, :) = xmean + sigma * randn(1, n) * diag(D) * B';
        else
            i = i + 1;
        end
    end
    
    % calculate fitness values
    fit = benchmark_func(x, func_num);
    
    % selection for offspring and parent
    [~, ind] = sort(fit);
    % find the index of best individual
    bestIdx = ind(1);
    bestFit = fit(bestIdx);  % the best fitness value
    bestX = x(bestIdx, :);  % the best individual
    
    if isfield(optionalArgs, 'fst_num')
        % get the best fitness value of the neighbour
        f_ij_b = bestFit;
        
        FEs = FEs + lambda;
        gens = gens + 1;
    else
        % update local variables
        FEs = FEs + lambda;
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
        
        xold = xmean;
        
        % recombination to generate new xmean
        xmean = weights * x(ind(1:mu), :);
        
        % update evolution paths: ps and pc
        ps = (1 - cs) * ps...
            + sqrt(cs * (2 - cs) * mueff) * (xmean - xold) * invsqrtC / sigma;
        %     hsig = sqrt(sum(ps.^2)) / sqrt(1-(1-cs)^(2 * FEs / lambda)) / chiN < 1.4 + 2/(n+1);
        hsig = sum(ps.^2) / (1-(1-cs)^(2 * FEs / lambda)) / n < 2 + 4/(n+1);
        pc = (1 - cc) * pc...
            + hsig * sqrt(cc * (2 - cc) * mueff)  * (xmean - xold) / sigma;
        
        % adapt covariance matrix C
        dif = (x(ind(1:mu), :) - repmat(xold, mu, 1)) / sigma;  % difference between x and xold
        C = (1 - c1 - cmu) * C + c1 * (pc' * pc) + cmu * dif' * diag(weights) * dif;
        
        % adapt step size sigma
        sigma = sigma * exp((norm(ps) / chiN - 1) * cs / damps);
        %     sigma = sigma * exp(min(1, (sqrt(sum(ps.^2)) / chiN - 1) * cs / damps));
        
        % update B and D from C
        C = triu(C) + triu(C, 1)';  % enforce symmetry
        [B, D] = eig(C);  % eigen decomposition, B == normalized eigenvectors
        D = sqrt(diag(D));  % D contains standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
        
        % build a struct array to save inner variables
        variables.FEs = FEs;
        variables.gens = gens;
        params = struct('B', B, 'D', D, 'xmean', xmean, 'sigma', sigma);
        variables.params = params; % save params in the first generation
        % save this generation to file
        saveEachGen(gens, saveGenPath, x, fit, variables);
        
        % ignore the last generation if: leftFEs < populationSize
        if options.MaxFEs - FEs < lambda
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

