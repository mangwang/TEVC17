function options = setOptions( varargin )
% Set initial parameters for all Evolutionary Algorithms.
%   Parameters:
%   varargin            - Variable-length input argument list
%                       [(key, value) argument pairs]
%   General Paramenters:
%   AlgoName            - The saving name
%                       [ string | 'GA' ]
%   PopulationType      - The type of Population being entered
%                       [ 'bitstring' | 'doubleVector' | 'bitstring' ]
%   PopInitRange        - Initial range of values a population may have
%                       [ row vector  | [-100, 100] ]
%   PopulationSize      - Positive scalar indicating the number of individuals
%                       [ positive scalar | 100 ]
%   EliteCount          - Number of best individuals that survive to next
%                         generation without any change
%                       [ positive scalar | 0.05*PopulationSize ]
%   CrossoverFraction   - The fraction of genes swapped between individuals
%                       [ positive scalar | 1 ]
%   Dim                 - Dimension of optimization functions
%                       [ positive scalar | 30 ]
%   MaxFEs              - Maximum number of Fitness Evaluations allowed
%                       [ positive scalar | 10000*Dim]
%   StallGenLimit    - Window size for termination criteria
%                       [ positive scalar | 30]
%   TolFun              - Termination tolerance on fitness value
%                       [ positive scalar | 1e-6 ]


% Init default options when no varargin given
if nargin == 0
    options = struct(...
        ... % general paramenters
        'AlgoName', 'GA', ...
        'PopulationType', 'doubleVector', ...
        'PopInitRange', [-100, 100], ...
        'PopulationSize', 100, ...
        'EliteCount', 0.05*100, ...
        'CrossoverFraction', 0.8, ...
        'Dim', 30, ...
        'MaxFEs', 1e4*30, ...
        'StallGenLimit', 50, ...
        'TolFun', 1e-6, ...
        ... % specific parameters
        'FitScalingType', 'rank', ...
        'ParentSelectionType', 'SUS', ...
        'TopSelectionQuantity', 0.4, ...
        'LinearSelectionMaximumSurvivalRate', 2, ...
        'TournamentSelectionTournamentSize', 4, ...
        'CrossoverType', 'wholeArithmetic', ...
        'CrossoverSBX_eta', 5, ...
        'MutationType', 'uniform', ...
        'MutationRate', 0.05, ...
        'MutationGaussianScale', 0.1, ...
        'MutationPolynomialEta', 20, ...
        'BoundConstraint', 'absorb' ...
        );
else
    % Set options, note that the first argument should be options struct
    options = varargin{1};
    for i = 2:2:(nargin-1)
        if isfield(options, varargin{i})
            options.(varargin{i}) = varargin{i+1};
        end
        if strcmp(varargin{i}, 'PopulationSize')
            options = setOptions(options, 'EliteCount', 0.05*options.PopulationSize);
        end
        if strcmp(varargin{i}, 'Dim')
            options = setOptions(options, 'MaxFEs', 1e4*options.Dim);
        end
    end
    if rem(i, 2)
        error('Variable-length input argument need to be odd!');
    end
    
end

end

