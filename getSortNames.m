function names = getSortNames(readPath)
% Read a folder, get all .mat file names sorted by integer.
%   Parameters:
%   readPath            - The path to read from
%                       [string]


f = dir([readPath, filesep, '*.mat']);
n = length(f);
count = zeros(1, n);
names = cell(1, n); % cell array to contain strings
for i = 1:n
    t = char(regexp(f(i).name, '\d+', 'match'));
    names{i} = t; % a string of integer
    count(i) = str2double(t);
end
[~, ind] = sort(count);
% reorder the names
names = names(ind);

