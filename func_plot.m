function fhead = func_plot(func_num)

% add benchmark path
addpath(genpath('benchmark_func_2013'));

if func_num == 14
    x = linspace(0, 2*pi, 100);
else
    x = linspace(-100, 100, 100);
end
y = x;

[X, Y] = meshgrid(x, y);
xx = reshape(X, length(x) ^ 2, 1);
yy = reshape(Y, length(y) ^ 2, 1);
f = benchmark_func([xx, yy], func_num);
f = reshape(f, length(x), length(x));

fig = figure();
fhead = surfc(Y,X,f);

print(fig, ['func_', num2str(func_num)], '-depsc');

end

