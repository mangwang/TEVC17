function [p_kw, H, Rmean, control_id, p_corrected, reject, best_ids] = ...
    kw_conover_control(groups, flag, alpha, method)
% Perfrom K-W Test followed by Connover post-hoc precedure and adjusted
% p-values procedure.
% Parameters:
% groups  - The input data
%           [matrix]  row is each measure, column is each group
%           [cell array]  each element of cell array is a group
% flag    - The flag to set the controlled group
%           [0 or 1]
%           if flag == 0, the minimum Rmean is set to the control;
%           if flag == 1, the maximum Rmean is set to the control;
% alpha   - significance level or FWER
%           [float]  0.05 is default
% method  - method for adjust p-values
%           [string]  'Holm' is default


% first perform the K-W test as the omnibus test
[p_kw, H, Rmean, ranked, N] = kw_test(groups);

if(~isvector(groups))
    k = size(groups, 2);
    sizes = NaN(1, k);
    for i = 1:k
        sizes(i) = size(groups, 1);
    end
else
    k = length(groups);
    sizes = NaN(1, k);
    for i = 1:k
        sizes(i) = length(groups{i});
    end
end

if(p_kw < alpha)
    k = length(sizes);
    s2 = (sum(ranked .^ 2) - N * (N+1) .^ 2 / 4) / (N - 1);
    if flag == 0
        [R_control, control_id] = min(Rmean);
    elseif flag == 1
        [R_control, control_id] = max(Rmean);
    end
    R_all = 1:k;
    R_all(control_id) = [];  % delete controlled group
    R_rest = R_all;
    
    % perform Conover post-hoc precedure
    R_diff = abs(Rmean(R_rest) - R_control) ./ sqrt(s2 * ((N-1-H) / (N-k)) .* ...
        (1 ./ sizes(R_rest) + 1 / sizes(control_id)));
    p_uncorrected = 1 - tcdf(R_diff, N-k);
    
    % adjusted p-values
    if(strcmp(method, 'Bonf'))
        p_corrected = p_uncorrected * (k-1);
        reject = (p_corrected <= alpha);
        best_ids = find(reject == 0);
        
    elseif(strcmp(method, 'Holm'))
        del_control = [1:(control_id-1) (control_id+1):k];
        [pc_uncorrected, ind] = sort(p_uncorrected, 'ascend');
        pc_corrected = pc_uncorrected .* (k - (1:(k-1)));
        p_corrected = zeros(size(pc_corrected));
        p_corrected(ind) = pc_corrected;
        reject = (p_corrected <= alpha);
        best_ids = [control_id del_control(reject == 0)];
        
    elseif(strcmp(method, 'Finn'))
        del_control = [1:(control_id-1) (control_id+1):k];
        [pc_uncorrected, ind] = sort(p_uncorrected, 'ascend');
        pc_corrected = 1 - (1-pc_uncorrected) .^ ((1:(k-1)) .\ (k-1));
        p_corrected = zeros(size(pc_corrected));
        p_corrected(ind) = pc_corrected;
        reject = (p_corrected <= alpha);
        best_ids = [control_id del_control(reject == 0)];
    end
    p_corrected(p_corrected > 1) = 1;
else
    control_id = 0;
    p_corrected = 1;
    reject = 0;
    best_ids = 1:length(sizes);
end
end
