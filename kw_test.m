function [p, H, Rmean, ranked, N] = kw_test(groups)
% Performs a non-parametric one-way
%   ANOVA to test the null hypothesis that independent samples from two or
%   more groups come from distributions with equal medians, and returns the
%   p-value for that test.
% If goups is a matrix, KRUSKALWALLIS treats each column as coming from a
%     separate group. This form of input is appropriate when each sample
%     has the same number of elements (balanced).
% If groups is a vector (or cell array), KRUSKALWALLIS treates each array
%     as coming from a separate group. This form of input is appropriate
%     when each sample has different number of elements (imbalanced).


if(~isvector(groups))
    [m, n] = size(groups);
    x = groups(:);
    N= length(x);
    [ranked, tieadj] = tiedrank(x);  % tie rank
    xr = reshape(ranked, m, n);
    Rm = sum(xr);
    Rmean = mean(xr);
    H  = 12 / (N * (N+1)) * sum(Rm .^ 2) / m - 3 * (N + 1);
    if (tieadj > 0)
        H = H / (1 - 2 * tieadj / (N .^ 3 - N));  % tie correction
    end
    df = n-1;
    p = 1 - chi2cdf(H, df);
else
    k = length(groups);
    
    x = groups{1};
    for i = 2:k
        x = [x groups{i}];
    end
    N = length(x);
    [ranked,tieadj] = tiedrank(x);  % tie rank
    
    inds = NaN(1,k+1);
    inds(1) = 0;
    Rm = NaN(1, k);
    Rmean = NaN(1,k);
    sizes = NaN(1, k);
    for i = 1:k
        isize = length(groups{i});
        sizes(i) = isize;
        inds(i+1)= inds(i) + isize;
        i_ind = (inds(i)+1):inds(i+1);
        irank = ranked(i_ind);
        Rm(i) = sum(irank);
        Rmean(i) = mean(irank);
    end
    
    H  = 12 / (N * (N+1)) * sum(Rm .^ 2 ./ sizes) - 3 * (N+1);
    if (tieadj > 0)
        H = H / (1 - 2 * tieadj / (N .^ 3 - N));  % tie correction
    end
    df = k-1;
    p = 1 - chi2cdf(H, df);
end