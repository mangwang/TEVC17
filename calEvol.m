function calEvol( run_num_normal, func_names, func_num, options )
% Calculate the features of population evolvability.
%   Parameters:
%   runNumNormal        - The number of normal runs
%                       [positive scalar]
%   func_names          - Function names
%                       [cell array of strings]
%   func_num            - The number of optimization functions
%                       [positive scalar]
%   options             - The options set by setOptions()
%                       [struct array]


nNeighbours = options.nNeighbours;  % the number of neighbours

% folder for save the best fitness value of all generated neighbours
saveNeiPath = ['result', filesep, 'neighBestFit', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
    'run_', num2str(run_num_normal)];

% folder of online sampling
readConvPath = ['result', filesep, 'conver_trend', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
    'run_', num2str(run_num_normal)];
convStr = load([readConvPath, filesep, 'bestFitSoFar.mat']);
FEsEachGen = convStr.FEsEachGen;

% how many geneartions for calculating indicators of evolvability
unifiedFEs = options.unifiedFEs;
nCalEvoPop = length(FEsEachGen(FEsEachGen <= unifiedFEs)) - 1;

if ~isdir(saveNeiPath)
    % folder for save the best fitness value of all generated neighbours
    saveImpPath = ['result', filesep, 'neighbours', filesep, 'dim_', num2str(options.Dim),...
        filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
        'run_', num2str(run_num_normal)];
    if ~isdir(saveImpPath)
        mkdir(saveImpPath);
    end
    optionalArgs.saveImpPath = saveImpPath; % save to optionalArgs
    
    % generate M neighbours (offspring) of each sampled population
    neighBestFit = zeros(nCalEvoPop, nNeighbours);  % init a matrix to save the best fitness value of all generated neighbours
    
    % path of sampled data
    readSamPath = ['result', filesep, 'raw_data', filesep, 'dim_', num2str(options.Dim),...
        filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
        'run_', num2str(run_num_normal)];
    
    parentFitStd  = zeros(nCalEvoPop, 1);
    for gen_num = 1:nCalEvoPop
        % load .mat data
        fileStr = load([readSamPath, filesep, num2str(gen_num), '.mat']);
        % load x and fit
        samPop = fileStr.x;
        samFit = fileStr.fit;
        % load variables
        optionalArgs.variables = fileStr.variables;
        % set other fields of optionalArgs
        optionalArgs.InitPopulation = samPop;
        optionalArgs.InitFitness = samFit;
        optionalArgs.fst_num = gen_num;
        optionalArgs.run_num_normal = run_num_normal;
        optionalArgs.nNeighbours = nNeighbours;
        
        % get std of parent fitness
        parentFitStd(gen_num) = std(samFit);
        
        algoName = options.AlgoName;
        for neigh_num = 1:nNeighbours
            % run the EA to generate one-step neighbours
            % the best fitness value of each generated neighbour is save to the 'saveImpPath' folder
            feval(algoName, neigh_num, func_names, func_num, options, optionalArgs);
            
            % load the best fitness value of each generated neighbour
            neighStr = load([saveImpPath, filesep, num2str(gen_num), '_', num2str(neigh_num), '.mat']);
            f_ij_b = neighStr.f_ij_b;
            % save to the matrix
            neighBestFit(gen_num, neigh_num) = f_ij_b;
        end
        
    end
    
    % load bestFitSoFar
    bestFitSoFar = convStr.bestFitSoFar;
    parentBestFitGlo = bestFitSoFar(1:nCalEvoPop)';  % load best fitness so far
    
    if ~isdir(saveNeiPath)
        mkdir(saveNeiPath);
    end
    save([saveNeiPath, filesep, 'neighBestFit.mat'], 'neighBestFit');
    save([saveNeiPath, filesep, 'parent.mat'], 'parentBestFitGlo', 'parentFitStd');
    
    % delete the temporal subfolder in the improve folder
    rmdir(saveImpPath, 's');
    
else
    % load parent and neighbours
    neighStr = load([saveNeiPath, filesep, 'neighBestFit.mat']);
    neighBestFit = neighStr.neighBestFit;
    parentStr = load([saveNeiPath, filesep, 'parent.mat']);
    parentBestFitGlo = parentStr.parentBestFitGlo;
    parentFitStd = parentStr.parentFitStd;
end

% ----- calculate indicators -----
neighBestFit = neighBestFit(1:nCalEvoPop, 1:nNeighbours);
% repmat values of parent
parentBestFitGlo = repmat(parentBestFitGlo(1:nCalEvoPop), 1, nNeighbours);
parentFitStd = repmat(parentFitStd(1:nCalEvoPop), 1, nNeighbours);

% calculate epp
impMat = neighBestFit < parentBestFitGlo;  % logical improvement matrix
impCntAllGens = sum(impMat, 2);  % imporoved count of all gens
eppAllGens = impCntAllGens ./ nNeighbours;  % epp of all gens

% calculate eap of each generation
impIdx = find(impCntAllGens >= 1);  % indices of gens that |N^+(P_i)| >= 1
% init eap
eapAllGens = zeros(nCalEvoPop, 1);  % init eapAllGens
for i = 1:length(impIdx)
    ig = impIdx(i);
    % calculate numerator
    numer = parentBestFitGlo(ig, impMat(ig, :)) - neighBestFit(ig, impMat(ig, :));
    % calculate denominator
    denom = parentFitStd(ig, impMat(ig, :));
%     denom(denom == 0) = eps;  % devide by zero correction
    denom(denom < 1e-8) = 1e-8;  % revise the denominator if it is a small value (threshold is set to be 1e-8)
    eapAllGens(ig) = mean(numer ./ denom);
end

% calculate eap and evp per FEs
FEsAllGens = FEsEachGen(2:(nCalEvoPop+1))' - FEsEachGen(1:nCalEvoPop)';
eapAllGens = eapAllGens ./ FEsAllGens;
evpAllGens = eppAllGens .* eapAllGens;

if ~isvector(evpAllGens)
    fprintf('The evpAllGens should be a vector!');
end
% save all the indicators
saveEvoPath = ['result', filesep, 'evolvability', filesep, 'dim_', num2str(options.Dim),...
    filesep, func_names{func_num}, filesep, options.AlgoName, filesep, ...
    'run_', num2str(run_num_normal)];
if ~isdir(saveEvoPath)
    mkdir(saveEvoPath);
end
save([saveEvoPath, filesep, 'indicators.mat'], 'eppAllGens', 'eapAllGens', 'evpAllGens');

end

