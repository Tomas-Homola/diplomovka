close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);


data = readmatrix("112x72.csv");
dataMono = readmatrix("mono_11x72_2018.csv");

%data(isnan(data)) = 0;
%dataMono(isnan(dataMono)) = 0;

dataAll = [data; dataMono];

n1 = 22;
n2 = 28;
n3 = 31;
n4 = 31;
n5 = 11;
N = n1 + n2 + n3 + n4 + n5;

C1start = 1;
C2start = n1 + 1;
C3start = n1 + n2 + 1;
C4start = n1 + n2 + n3 + 1;
C5start = n1 + n2 + n3 + n4 + 1;

C1end = n1;
C2end = n1 + n2;
C3end = n1 + n2 + n3;
C4end = n1 + n2 + n3 + n4;
C5end = n1 + n2 + n3 + n4 + n5;

ind = [C1start, C1end; C2start, C2end; C3start, C3end; C4start, C4end; C5start, C5end];

D = pdist(dataAll,"jaccard");
D = squareform(D);

Y = cmdscale(D);

figure
hold on
plot(Y(ind(1,1):ind(1,2), 1), Y(ind(1,1):ind(1,2), 2), 'r.', 'MarkerSize', 15);
plot(Y(ind(2,1):ind(2,2), 1), Y(ind(2,1):ind(2,2), 2), 'g.', 'MarkerSize', 15);
plot(Y(ind(3,1):ind(3,2), 1), Y(ind(3,1):ind(3,2), 2), 'b.', 'MarkerSize', 15);
plot(Y(ind(4,1):ind(4,2), 1), Y(ind(4,1):ind(4,2), 2), 'm.', 'MarkerSize', 15);
plot(Y(ind(5,1):ind(5,2), 1), Y(ind(5,1):ind(5,2), 2), 'k.', 'MarkerSize', 15);
hold off

%%
temp = data(18,:)';
dist = norm(temp - temp);

%%
for i = 1:112
	fprintf("%d: %d\n", i, nnz(isnan(data(i,:))));
end




