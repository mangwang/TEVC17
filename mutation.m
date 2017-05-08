function mutationKids = mutation(parents, options)
% Perform mutation on the parents.
% Each gene of one child has mutationRate chance to be mutated,
% and at least one gene need to be mutated.
%   Parameters:
%   parents             - genotype of the selected parents
%                       [Matrix of Dim*nParents dimension]
%   options             - options
%                       [struct array]


mutationKids = zeros(size(parents));
mutationRate = options.MutationRate;

range = options.PopInitRange;
lb = range(1);
ub = range(2);

if(strcmpi(options.PopulationType, 'doubleVector'))
    switch options.MutationType
        case 'uniform'
            % Reference:
            % Agoston E. Eiben and J. E. Smith. 2015.
            % Introduction to Evolutionary Computing, 2nd Edition. SpringerVerlag.
            % Page 57.
            
            for i = 1:size(parents, 1)
                child = parents(i, :);
                % Each gene of one child has mutationRate chance to be mutated
                mutationPoints = find(rand(1, length(child)) < mutationRate);
                %                 % At least one one gene need to be mutated
                %                 if(isempty(mutationPoints))
                %                     mutationPoints = ceil(rand * length(child));
                %                 end
                % each gene is replaced with a value chosen randomly from the range
                child(mutationPoints) = lb + rand(size(mutationPoints)) .* (ub - lb);
                
                mutationKids(i, :) = child;
            end
            
        case 'gaussian'
            % Reference:
            % Agoston E. Eiben and J. E. Smith. 2015.
            % Introduction to Evolutionary Computing, 2nd Edition. SpringerVerlag.
            % Page 57.
            
            for i = 1:size(parents, 1)
                child = parents(i, :);
                % Each gene of one child has mutationRate chance to be mutated
                mutationPoints = find(rand(1, length(child)) < mutationRate);
                %                 % At least one one gene need to be mutated
                %                 if(isempty(mutationPoints))
                %                     mutationPoints = ceil(rand * length(child));
                %                 end
                % each gene is replaced with a value chosen randomly from the range
                sigma = options.MutationGaussianScale .* (ub - lb);
                child(mutationPoints) = child(mutationPoints) + ...
                    sigma .* randn(size(mutationPoints));
                
                mutationKids(i,:) = child;
            end
        case 'polynomial'
            % Reference:
            % Deb, Kalyanmoy, and Debayan Deb. "Analysing mutation schemes
            % for real-parameter genetic algorithms." International Journal
            % of Artificial Intelligence and Soft Computing 4.1 (2014): 1-28.
            
            for i = 1:size(parents, 1)
                child = parents(i, :);
                % Each gene of one child has mutationRate chance to be mutated
                mutationPoints = find(rand(1, length(child)) < mutationRate);
                %                 % At least one gene need to be mutated
                %                 if(isempty(mutationPoints))
                %                     mutationPoints = ceil(rand * length(child));
                %                 end
                % each gene is replaced with a value chosen randomly from the range
                eta_m = options.MutationPolynomialEta;
                
                delta = zeros(size(mutationPoints));
                u = rand(size(mutationPoints));
                delta(u <= 0.5) = (2 .* u(u <= 0.5)) .^ (1 / (1 + eta_m)) - 1;
                delta(u > 0.5) = 1 - (2 .* (1 - u(u > 0.5))) .^ (1 / (eta_m + 1));
                
                child(mutationPoints(u <= 0.5)) = child(mutationPoints(u <= 0.5)) ...
                    + delta(u <= 0.5) .* (child(mutationPoints(u <= 0.5)) - lb);
                child(mutationPoints(u > 0.5)) = child(mutationPoints(u > 0.5)) ...
                    + delta(u > 0.5) .* (ub - child(mutationPoints(u > 0.5)));
                
                mutationKids(i, :) = child;
            end
        otherwise
            error('Wrong Mutation Type!')
    end
elseif(strcmpi(options.PopulationType, 'bitstring'))
    error('Need to be finished!');
else
    error('Wrong Population Type!')
end


end