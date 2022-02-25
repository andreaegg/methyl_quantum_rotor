function [] = xyz_file2plot(filename,path,headerlines)

all = 1e10;

if nargin < 2
    path = '';
    headerlines = 0;
end

if nargin < 3 
    headerlines = 0;
end

input = strcat(path,filename);
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

figure(100),clf
for m = 1:length(bond)
    plot3([xyz(connect(m,1),1) xyz(connect(m,2),1)],[xyz(connect(m,1),2) xyz(connect(m,2),2)],[xyz(connect(m,1),3) xyz(connect(m,2),3)],'k-');
    hold on;
    box on; grid on;
end
for m = 1:length(number)
    plot3(xyz(m,1),xyz(m,2),xyz(m,3),'.','Color',atomcolors(atom{m}),'MarkerSize',25)
    text(xyz(m,1)+0.1,xyz(m,2)+0.1,xyz(m,3)+0.1,num2str(number(m)));
end
xlabel('X coord')
ylabel('Y coord')
zlabel('Z coord')

end
