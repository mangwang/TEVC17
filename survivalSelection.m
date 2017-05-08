function parents = survivalSelection( fitness, nParents, surSelType )
% Perform the survival selection to generate the offspring.
% The returned parents is the indices of chosen individuals.
%   Parameters:
%   fitness             - fitness values
%                       [column vector]
%   nParents            - number of parents for next generation
%                       [positive scalar]
%   surSelType          - survical selection type
%                       [string]


fitness = fitness(:); % fitness values of current population
nPop = length(fitness);

switch surSelType
    case 'comma'
        % Perform (mu, lambda) selection
        [~, ind] = sort(fitness);
        parents = ind(1:nParents);
    case 'roundRobin'
        % set tournament size q
        q = 10;
        % select competitors
        selCompes = fitness(randperm(nPop, q));
        % start competition
        score = zeros(nPop, 1);
        for i = 1:nPop
            score(i) = sum(fitness(i) < selCompes);
        end
        % select best mu individuals
        [~, ind] = sort(score, 'descend');
        parents = ind(1:nParents);
    otherwise
        error('Wrong Parent Selection Type!');
end

end

