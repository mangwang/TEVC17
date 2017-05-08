function [ lb, ub ] = get_lb_ub( func_num )
% Get lower bound and upper bound of each test function.
%   Parameters:
%   func_num            - Number of the test function
%                       [positive scalar]


% get lb and ub of the whole search space
if func_num >= 1 && func_num <= 13
    lb = -100;
    ub = 100;
end

if func_num == 14
    lb = 0;
    ub = 2*pi;
end

end

