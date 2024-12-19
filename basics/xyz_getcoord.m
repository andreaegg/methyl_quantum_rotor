function [xyz,atom,connect] = xyz_getcoord(filename,path,headerlines)

all = 1e10;

if nargin < 2
    path = '';
    headerlines = 0;
end

if nargin < 3 
    headerlines = 0;
end

if isempty(path)
    input = filename;
else
    input = strcat(path,'\',filename);
end

data = importdata(input,'\t',all);

% remove headerlines
data = data(1+headerlines:end);
data = strtrim(data);
xyz  = zeros(numel(data),3);

for k = 1:numel(data)
    curratom = split(convertCharsToStrings(data{k}));
    atom(k)  = curratom(1);
    xyz(k,1) = str2num(curratom(2));
    xyz(k,2) = str2num(curratom(3));
    xyz(k,3) = str2num(curratom(4));
end

number = 1:1:numel(data);
[connect,bond] = atomsconnectivity(atom);

end
