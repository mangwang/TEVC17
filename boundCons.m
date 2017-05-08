function x = boundCons( x, options )
% Boundary constraint handling scheme.
%   Parameters:
%   x                   - The genotype of each generation
%                       [Matrix]
%   options             - The options set by setOptions()
%                       [struct array]


% get lower and upper bound
lb = options.PopInitRange(1);
ub = options.PopInitRange(2);

% update bound using one of the bound constraint scheme
switch options.BoundConstraint
    % Reference:
    % Gandomi, Amir Hossein, and Xin-She Yang. 
    % "Evolutionary boundary constraint handling scheme." 
    % Neural Computing and Applications 21.6 (2012): 1449-1462..
    
    case 'absorb'
        x = max(min(x, ub), lb);
    case 'random'
        idx = x > ub | x < lb;
        x(idx) = lb + rand(length(x(idx)), 1) .* (ub - lb);
    case 'reflect'
        % outside of the lower bound
        idxOutLower = x < lb;
        x(idxOutLower) = lb - (x(idxOutLower) - lb);
        x(idxOutLower) = min(x(idxOutLower), ub);
        % outside of the upper bound
        idxOutUpper = x > ub;
        x(idxOutUpper) = ub - (x(idxOutUpper) - ub);
        x(idxOutUpper) = max(x(idxOutUpper), lb);
    case 'toroidal'
        % outside of the lower bound
        idxOutLower = x < lb;
        x(idxOutLower) = ub + (x(idxOutLower) - lb);
        x(idxOutLower) = max(x(idxOutLower), lb);
        % outside of the upper bound
        idxOutUpper = x > ub;
        x(idxOutUpper) = lb + (x(idxOutUpper) - ub);
        x(idxOutUpper) = min(x(idxOutUpper), ub);
end


end

