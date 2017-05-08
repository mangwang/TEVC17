% Execute the whole process
clear; clc; close all;

% add benchmark path
addpath(genpath('benchmark_func_2013'));

% set parameters
func_names = {'1-Sphere', '2-Ellipsoidal', '3-BentCigar', '4-Discus', '5-DiffPowers', ...
    '6-Rosenbrock', '7-SchafferF7', '8-Ackley', '9-Weierstrass', '10-Griewank', ...
    '11-Rastrigin', '12-Schwefel', '13-Katsuura', '14-Radar'
    };

algo_names = {'GA', 'CMA_ES', 'CoDE', 'SPSO2011', 'ABC'};
dims = [2, 10, 30, 50];

runStart = 1;
runEnd = 1;


% --------------------  Visualization of Generational Measures with Convergence Curve -------------------- %

% Online sampling using each algorithm until reaching the limFEs
options = setOptions(); % init options
optionalArgs = struct();
for d = 3
    dim = dims(d);
    for func_num = 1:14
        [lb, ub] = get_lb_ub(func_num);
        
        if func_num == 14
            dim = 20;
        end
        
        for algo_num = 1:5
            algo_name = algo_names{algo_num};
            
            % set population size (default is 100) and other params
            if strcmp(algo_name, 'CMA_ES')
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 4+floor(3*log(dim)), 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'CoDE')
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 30, 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'SPSO2011')
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 40, 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'ABC')
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 50, 'Dim', dim, 'PopInitRange', [lb, ub]);
            else
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 100, 'Dim', dim, 'PopInitRange', [lb, ub]);
            end
            
            % set the limFEs to maxFEs/5
            options = setOptions(options, 'MaxFEs', options.MaxFEs / 5);
            
            % run an EA for online sampling
            for run_num = runStart:runEnd
                feval(algo_name, run_num, func_names, func_num, options, optionalArgs);
            end
            
        end
    end
end


% Calculate evolvability
options = setOptions(); % init options
for d = 3
    dim = dims(d);
    
    for func_num = 1:14
        [lb, ub] = get_lb_ub(func_num);
        
        if func_num == 14
            dim = 20;
        end
        
        % calculate evolvability
        for algo_num = 1:5
            algo_name = algo_names{algo_num};
            
            % set population size (default is 100) and other params
            if strcmp(algo_name, 'CMA_ES')
                options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 4+floor(3*log(dim)), 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'CoDE')
                options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 30, 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'SPSO2011')
                options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 40, 'Dim', dim, 'PopInitRange', [lb, ub]);
            elseif strcmp(algo_name, 'ABC')
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 50, 'Dim', dim, 'PopInitRange', [lb, ub]);
            else
                options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 100, 'Dim', dim, 'PopInitRange', [lb, ub]);
            end
            options.unifiedFEs = options.MaxFEs / 5;
            options.nNeighbours = 5;
            
            for run_num = runStart:runEnd
                calEvol(run_num, func_names, func_num, options);
            end
            
        end
    end
end


% Visulization
options = setOptions(); % init options
for d = 3
    dim = dims(d);
    for func_num = 1:14
        [lb, ub] = get_lb_ub(func_num);
        
        if func_num == 14
            dim = 20;
        end
        
        for run_num = runStart:runEnd
            
            for algo_num = 1:5
                algo_name = algo_names{algo_num};
                
                % set population size (default is 100) and other params
                if strcmp(algo_name, 'CMA_ES')
                    options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 4+floor(3*log(dim)), 'Dim', dim, 'PopInitRange', [lb, ub]);
                elseif strcmp(algo_name, 'CoDE')
                    options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 30, 'Dim', dim, 'PopInitRange', [lb, ub]);
                elseif strcmp(algo_name, 'SPSO2011')
                    options = setOptions(options, 'AlgoName', algo_name,  'PopulationSize', 40, 'Dim', dim, 'PopInitRange', [lb, ub]);
                elseif strcmp(algo_name, 'ABC')
                    options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 50, 'Dim', dim, 'PopInitRange', [lb, ub]);
                else
                    options = setOptions(options, 'AlgoName', algo_name, 'PopulationSize', 100, 'Dim', dim, 'PopInitRange', [lb, ub]);
                end
                
                % plot
                options.form = 'eps';  % saved format
                visulization(run_num, func_names, func_num, options);
            end
            
        end
    end
end



% --------------------  Algorithm Selection Framework  -------------------- %

% modify the unified FEs if an algorithm achieves the best performance among
% all candidate algorithms but used less FEs
uniFEsEachFunc = zeros(14, 1);
for d = 3
    dim = dims(d);
    for func_num = 1:14
        
        if func_num == 14
            dim = 20;
        end
        
        % find the best fitness value achieved on each function, and then unified the FEs
        % of all algorithms to the minimal FEs that the algorithm first
        % reached this best fitness value (in float-point precision)
        bestFitEachFunc = inf;
        unifiedFEs = inf;
        for algo_num = 1:5
            algo_name = algo_names{algo_num};
            saveAlgoName = strrep(algo_name, '_', '-');  % replace _ to - for saving into file
            
            % read the global best fitness values of online sampling
            for run_num_normal = runStart:runEnd
                readConvPath = ['result', filesep, 'conver_trend', filesep, 'dim_', num2str(dim),...
                    filesep, func_names{func_num}, filesep, algo_name, filesep, ...
                    'run_', num2str(run_num_normal)];
                convStr = load([readConvPath, filesep, 'bestFitSoFar.mat']);
                bestFitSoFar = convStr.bestFitSoFar;
                FEsEachGen = convStr.FEsEachGen;
                % performing the float-point correction (let (number - floor(number) < 1e-8) = 0)
                maskbestFit = (bestFitSoFar - floor(bestFitSoFar) < 1e-8 & bestFitSoFar - floor(bestFitSoFar) > 0);
                if any(maskbestFit)
                    % do the correction
                    bestFitSoFar(maskbestFit) = floor(bestFitSoFar(maskbestFit));
                end
                % find the FEs that the algorithm first reached this fitness
                bestFitLast = bestFitSoFar(end);
                idx = find(bestFitSoFar == bestFitLast, 1);
                usedFEs = FEsEachGen(idx);
                
                % set the unified FEs
                if bestFitLast < bestFitEachFunc  % best fitness value
                    bestFitEachFunc = bestFitLast;
                    unifiedFEs = usedFEs;
                    bestAlgo = saveAlgoName;
                elseif bestFitLast == bestFitEachFunc
                    if usedFEs < unifiedFEs  % fastest
                        unifiedFEs = usedFEs;
                        bestAlgo = saveAlgoName;
                    end
                end
                
            end
        end
        uniFEsEachFunc(func_num) = unifiedFEs;
    end
end
% save the uniFEsEachFunc to mat file
saveuniFEsPath = ['result', filesep, 'FEs_stats'];
if ~isdir(saveuniFEsPath)
    mkdir(saveuniFEsPath);
end
save([saveuniFEsPath, filesep, 'uniFEsEachFunc.mat'], 'uniFEsEachFunc');


% Get the data of evolvability
for d = 3
    dim = dims(d);
    for func_num = 1:14
        
        if func_num == 14
            dim = 20;
        end
        
        saveEvpPath = ['result', filesep, 'stats', filesep, 'data', filesep, 'evp', filesep,...
            'dim_', num2str(dim), filesep, 'f_', num2str(func_num)];
        if ~isdir(saveEvpPath)
            mkdir(saveEvpPath);
        end
        
        for run_num_normal = runStart:runEnd
            
            evpAllAlgos = cell(1);  % init evp of all algorithms
            for algo_num = 1:5
                algo_name = algo_names{algo_num};
                
                loadEvoPath = ['result', filesep, 'evolvability', filesep, 'dim_', num2str(dim),...
                    filesep, func_names{func_num}, filesep, algo_name, filesep, ...
                    'run_', num2str(run_num_normal)];
                indStr = load([loadEvoPath, filesep, 'indicators.mat']);
                evpAllGens = indStr.evpAllGens;
                
                if exist('uniFEsEachFunc', 'var')
                    unifiedFEs = uniFEsEachFunc(func_num);  % set unified FEs
                else
                    uniFEsStr = load(['result', filesep, 'FEs_stats', filesep, 'uniFEsEachFunc.mat']);
                    unifiedFEs = uniFEsStr.uniFEsEachFunc(func_num);  % set unified FEs
                end
                
                % folder of online sampling
                readConvPath = ['result', filesep, 'conver_trend', filesep, 'dim_', num2str(dim),...
                    filesep, func_names{func_num}, filesep, algo_name, filesep, ...
                    'run_', num2str(run_num_normal)];
                convStr = load([readConvPath, filesep, 'bestFitSoFar.mat']);
                FEsEachGen = convStr.FEsEachGen;
                
                % real calculated number of generations for estimating population evolvability
                nCalEvoPop = length(FEsEachGen(FEsEachGen <= unifiedFEs)) - 1;
                
                evpAllAlgos{algo_num} = evpAllGens(1:nCalEvoPop)';
            end
            % save the cell array of data
            save([saveEvpPath, filesep, 'evp.mat'], 'evpAllAlgos');
        end
    end
end


% Get the data of benchmark
for d = 3
    dim = dims(d);
    for func_num = 1:14
        
        if func_num == 14
            dim = 20;
        end
        
        saveBenchPath = ['result', filesep, 'stats', filesep, 'data', filesep, 'bench', filesep,...
            'dim_', num2str(dim), filesep, 'f_', num2str(func_num)];
        if ~isdir(saveBenchPath)
            mkdir(saveBenchPath);
        end
        
        benchAllAlgos = cell(1);  % init performance of all aogorithms
        for algo_num = 1:5
            algo_name = algo_names{algo_num};
            loadBenchBase = ['result', filesep, 'benchmark', filesep, filesep, 'conver_trend_all_runs', ...
                filesep, 'dim_', num2str(dim), filesep, func_names{func_num}, filesep, algo_name];
            benchStr = load([loadBenchBase, filesep, 'bestFitAllRuns.mat']);
            bestFitAllRuns = benchStr.bestFitAllRuns;
            finalPerform = bestFitAllRuns(:, end);  % final perfermance of all runs
            benchAllAlgos{algo_num} = finalPerform';
            % save the cell array of data
            save([saveBenchPath, filesep, 'bench.mat'], 'benchAllAlgos');
            
        end
    end
end


% Statistical analysis
% open a file for writing results of evp
evp_res = fopen(['result', filesep, 'evp.csv'], 'w');
% open a file for writing results of benchmark
evp_bench = fopen(['result', filesep, 'bench.csv'], 'w');

for d = 3
    dim = dims(d);
    for func_num = 1:14
        
        if func_num == 14
            dim = 20;
        end
        
        % load evp groups
        loadEvpPath = ['result', filesep, 'stats', filesep, 'data', filesep, 'evp', filesep,...
            'dim_', num2str(dim), filesep, 'f_', num2str(func_num)];
        evpAllAlgosStr = load([loadEvpPath, filesep, 'evp.mat']);
        evpAllAlgos = evpAllAlgosStr.evpAllAlgos;
        
        % perform statistical test of evp sequences
        [~, ~, Rmean, ~, ~, ~, best_ids] = ...
            kw_conover_control(evpAllAlgos, 1, 0.05, 'Holm');
        
        k = length(Rmean);  % number of candidate algorithms
        for i = 1:k
            if ismember(i, best_ids)
                if i ~= k
                    fprintf(evp_res, '%d (*),', round(Rmean(i)));
                else
                    fprintf(evp_res, '%d (*)\n', round(Rmean(i)));
                end
            else
                if i ~= k
                    fprintf(evp_res, '%d,', round(Rmean(i)));
                else
                    fprintf(evp_res, '%d\n', round(Rmean(i)));
                end
            end
        end
        
        % load benchmark groups
        loadBenchPath = ['result', filesep, 'stats', filesep, 'data', filesep, 'bench', filesep,...
            'dim_', num2str(dim), filesep, 'f_', num2str(func_num)];
        benchAllAlgosStr = load([loadBenchPath, filesep, 'bench.mat']);
        benchAllAlgos = benchAllAlgosStr.benchAllAlgos;
        
        % perform statistical test of benchmark sequences
        [~, ~, Rmean, ~, ~, ~, best_ids] = ...
            kw_conover_control(benchAllAlgos, 0, 0.05, 'Holm');
        
        k = length(Rmean);  % number of candidate algorithms
        for i = 1:k
            if ismember(i, best_ids)
                if i ~= k
                    fprintf(evp_bench, '%d (*),', round(Rmean(i)));
                else
                    fprintf(evp_bench, '%d (*)\n', round(Rmean(i)));
                end
            else
                if i ~= k
                    fprintf(evp_bench, '%d,', round(Rmean(i)));
                else
                    fprintf(evp_bench, '%d\n', round(Rmean(i)));
                end
            end
        end
        
    end
end
fclose(evp_res);
fclose(evp_bench);

