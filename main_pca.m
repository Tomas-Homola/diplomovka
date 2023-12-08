clear
clc

%% Data import and prep
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

% INFO: nacitanie dat, najprv Sentinel data, potom point cloud data
% sentinel data
dataSent = readmatrix("112x72.csv");
dataMonoSent = readmatrix("mono_11x72_2018.csv");

% point cloud data
% dataRM = readmatrix("RM_curves_rev2.3.csv");
dataRM = readmatrix("RM_rev2.3_onlyVegHigh.csv");
dataMonoRM = readmatrix("RM_mono_onlyVegHigh.csv");

% INFO: tu sa vybera, ktore metriky z point cloudu zahrnut
% 18 19 20 21
index = @(i) ((i-1)*4 + 1):1:(i*4);

selected = [];
% INFO: ked chceme vsetky, tak treba dat [ 1:1:22 ], alebo len konkretne cisla
% INFO: ocislovane metriky su v subore dataFeatureExtractorCurve.m riadky 123-146
metrics = [5 13 20]; 
for i = metrics
	selected = [ selected index(i) ];
end

% INFO: vybrat moznost podla toho, ci chceme iba Sentinel, iba point cloud alebo oboje
data = dataSent;
% data = dataRM(:, selected);
% data = [dataSent dataRM(:, selected)];

dataMono = dataMonoSent;
% dataMono = dataMonoRM(:, selected);
% dataMono = [dataMonoSent dataMonoRM(:, selected)];

% INFO: normovanie dat pred PCA
maxColumnValuesBeforePCA = max(data);
minColumnValuesBeforePCA = min(data);
maxColumnValuesBeforePCA(maxColumnValuesBeforePCA == 0) = 1;

dataNormed = (data - minColumnValuesBeforePCA) ./ (maxColumnValuesBeforePCA - minColumnValuesBeforePCA);
dataNormed(:,maxColumnValuesBeforePCA == minColumnValuesBeforePCA) = 0;

% INFO: tu sa spravi PCA
[coef, score, ~, ~, explained, ~] = pca(dataNormed, "Algorithm","eig");
fprintf("explained: %.2f\n", sum(explained(1:3)));

% INFO: potom sa znormuju nove suradnice (premenna "score")
% INFO: preco sa to normuje aj po PCA si uz nepamatam, ale tak vzdy sa da
% tato cast kodu zakomentovat ...
maxColumnValuesAfterPCA = max(score);
maxColumnValuesAfterPCA(maxColumnValuesAfterPCA == 0) = 1;
minColumnValuesAfterPCA = min(score);

score = (score - minColumnValuesAfterPCA) ./ (maxColumnValuesAfterPCA - minColumnValuesAfterPCA);
score(isnan(score)) = 0;

% INFO: normovanie dat pre monokultury podla min/max z biotopov
dataMonoNormed = (dataMono - minColumnValuesBeforePCA) ./ (maxColumnValuesBeforePCA - minColumnValuesBeforePCA);
dataMonoNormed(:, maxColumnValuesBeforePCA == minColumnValuesBeforePCA) = 0;

% INFO: PCA pre monokultury
scoreMono = dataMonoNormed * coef;

% INFO: normovanie monokultur po PCA
% INFO: ... a potom aj toto zakomentovat, ak sa po PCA nechce normovat
scoreMono = (scoreMono - minColumnValuesAfterPCA) ./ (maxColumnValuesAfterPCA - minColumnValuesAfterPCA);
% scoreMono(isnan(scoreMono)) = 0;

% INFO: pomocne indexi
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

z = 3;
X = [score(:,1); scoreMono(:,1)];
Y = [score(:,2); scoreMono(:,2)];
Z = [score(:,z); scoreMono(:,z)];

figure
hold on
scatter(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), 20, "red", "filled");
scatter(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), 20, "green", "filled");
scatter(X(ind(3,1):ind(3,2)), Y(ind(3,1):ind(3,2)), 20, "blue", "filled");
scatter(X(ind(4,1):ind(4,2)), Y(ind(4,1):ind(4,2)), 20, "magenta", "filled");
scatter(X(ind(5,1):ind(5,2)), Y(ind(5,1):ind(5,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')

figure
scatter3(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), Z(ind(1,1):ind(1,2)), 20, "red", "filled");
hold on
scatter3(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), Z(ind(2,1):ind(2,2)), 20, "green", "filled");
scatter3(X(ind(3,1):ind(3,2)), Y(ind(3,1):ind(3,2)), Z(ind(3,1):ind(3,2)), 20, "blue", "filled");
scatter3(X(ind(4,1):ind(4,2)), Y(ind(4,1):ind(4,2)), Z(ind(4,1):ind(4,2)), 20, "magenta", "filled");
scatter3(X(ind(5,1):ind(5,2)), Y(ind(5,1):ind(5,2)), Z(ind(5,1):ind(5,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')

% Parallel coordinates
% INFO: vypocet priemerov a std
means = cell(5,1);
sigmas = cell(5,2);

for i = 1:4
	means{i} = mean( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,1} = means{i} - std( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,2} = means{i} + std( score( ind(i,1):ind(i,2), : ) );
end

means{5} = mean( scoreMono( 1:11, : ) );
sigmas{5,1} = means{5} - std( scoreMono( 1:11, : ) );
sigmas{5,2} = means{5} + std( scoreMono( 1:11, : ) );

% INFO: vykreslenie
figure; 
hold on

features = size(score,2);

% INFO: toto je tu len kvoli legende
plot([0,0],[0,0],'-r');
plot([0,0],[0,0],'-g');
plot([0,0],[0,0],'-b');
plot([0,0],[0,0],'-m');
plot([0,0],[0,0],'-k');

% 91E0
plot(1:1:features,  means{1}(1:features),  '.-r', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{1,1}(1:features),':r',  "LineWidth", 2)
plot(1:1:features, sigmas{1,2}(1:features),':r',  "LineWidth", 2)

% 91F0
plot(1:1:features,  means{2}(1:features),  '.-g', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{2,1}(1:features),':g',  "LineWidth", 2)
plot(1:1:features, sigmas{2,2}(1:features),':g',  "LineWidth", 2)

% 91G0
plot(1:1:features,  means{3}(1:features),  '.-b', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{3,1}(1:features),':b',  "LineWidth", 2)
plot(1:1:features, sigmas{3,2}(1:features),':b',  "LineWidth", 2)

% 9110
plot(1:1:features,  means{4}(1:features),  '.-m', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{4,1}(1:features),':m',  "LineWidth", 2)
plot(1:1:features, sigmas{4,2}(1:features),':m',  "LineWidth", 2)

% monokultury
plot(1:1:features,  means{5}(1:features),  '.-k', "LineWidth", 3, "MarkerSize", 20)
plot(1:1:features, sigmas{5,1}(1:features),':k',  "LineWidth", 1.5)
plot(1:1:features, sigmas{5,2}(1:features),':k',  "LineWidth", 1.5)

legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
hold off

% INFO: povodne vykreslenie vsetkych bodov
% figure; 
% hold on
% % features = length(selected);
% features = length(coef);
% 
% for i = C1start:C1end 
% 	plot(1:1:features, score(i,1:features), '.-r');
% end
% for i = C2start:C2end 
% 	plot(1:1:features, score(i,1:features), '.-g');
% end
% for i = C3start:C3end
% 	plot(1:1:features, score(i,1:features), '.-b');
% end
% for i = C4start:C4end 
% 	plot(1:1:features, score(i,1:features), '.-m');
% end
% for i = 1:11
% 	plot(1:1:features, scoreMono(i,1:features), '.-k'); 
% end
% % legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
% hold off



