function mutMatrix2ARM(file_path, pname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
input_txt = strcat(file_path, pname, '_mutMatrix.txt');
data = csvread(input_txt);

data_2 = data; % if data is not bool, transform to bool (data_2)
% data_2(data_2<=5)=0;
% data_2(data_2>5)=1;
% csvwrite('clrmarix2002_bool.csv',data_2);

%--- transform to "transaction data" for ARM (association rule mining)
L = size(data_2,1);
Indx = cell(L,1);
for i=1:L
    row = data_2(i,:);
    ind = find(row==1);
    Indx{i} = ind;
end
%--- write to file:
filename = strcat(file_path, pname, '_forARM', '.csv');
cell2csv(filename,Indx,',')
end

