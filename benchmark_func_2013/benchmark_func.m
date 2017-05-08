function fit = benchmark_func(x, func_num)

if func_num >=1 && func_num <=13
    % call cec13 function
    fit = cec13_func(x', func_num);
    fit = fit';
end

if func_num == 14
    % spread spectrum radar polyphase problem
    % D = 20, bound=[0,2*pi]
    
    [ps, d] = size(x);
    fit = zeros(ps, 1);
    
    for p = 1:ps
        m = 2 * d - 1;
        M = 2 * m;
        hsum = zeros(1, M);
        for kk = 1:m
            if rem(kk,2) ~= 0
                i= (kk+1) / 2;
                hsum(kk) = 0;
                for j = i:d
                    summ = sum(x(p, (abs(2*i-j-1)+1):j));
                    hsum(kk)= hsum(kk) + cos(summ);
                end
            else
                i = kk / 2;
                hsum(kk) = 0;
                for j = (i+1):d
                    summ = sum(x(p, (abs(2*i-j)+1):j));
                    hsum(kk)= hsum(kk) + cos(summ);
                end
                hsum(kk) = hsum(kk) + 0.5;
            end
        end
        hsum(m+1:M) = -hsum(1:m);
        fit(p) = max(hsum);
    end
    
end

end

