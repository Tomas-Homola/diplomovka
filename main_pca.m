%% Data import and prep
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

dataSent = readmatrix("112x72.csv");
dataRM = readmatrix("RM_curves_rev2.3.csv");
dataMono = readmatrix("RM_curvesMono_inSVK.csv");

% 18 19 20 21
selected = [index(18) index(19) index(20) index(21) index(22)];

data = [dataSent dataRM(:, selected)];

maxColumnValuesBeforePCA = max(data);
minColumnValuesBeforePCA = min(data);
maxColumnValuesBeforePCA(maxColumnValuesBeforePCA == 0) = 1;

dataNormed = (data - minColumnValuesBeforePCA) ./ (maxColumnValuesBeforePCA - minColumnValuesBeforePCA);
dataNormed(:,maxColumnValuesBeforePCA == minColumnValuesBeforePCA) = 0;

% dataNormed2 = normalize(data, "range");

% PCA
[coef, score, ~, ~, ~, ~] = pca(dataNormed, "Algorithm","eig");

maxColumnValuesAfterPCA = max(score);
maxColumnValuesAfterPCA(maxColumnValuesAfterPCA == 0) = 1;
minColumnValuesAfterPCA = min(score);

score = (score - minColumnValuesAfterPCA) ./ (maxColumnValuesAfterPCA - minColumnValuesAfterPCA);
score(:,maxColumnValuesAfterPCA == minColumnValuesAfterPCA) = 0;


% dataMonoNormed = (dataMono - minColumnValuesBeforePCA) ./ (maxColumnValuesBeforePCA - minColumnValuesBeforePCA);
% dataMonoNormed(:, maxColumnValuesBeforePCA == minColumnValuesBeforePCA) = 0;
% 
% % centrovanie novych bodov
% dataMonoNormed = dataMonoNormed - mean(dataMonoNormed);
% 
% % PCA pre nove body
% scoreMono = dataMonoNormed * coef;
% 
% % normovanie po PCA
% scoreMono = (scoreMono - minColumnValuesAfterPCA) ./ (maxColumnValuesAfterPCA - minColumnValuesAfterPCA);


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



% X = [score(:,1); scoreMono(:,1)];
% Y = [score(:,2); scoreMono(:,2)];

X = score(:,1);
Y = score(:,2);

figure
scatter(X(C1start:C1end), Y(C1start:C1end), 20, "red", "filled");
hold on
scatter(X(C2start:C2end), Y(C2start:C2end), 20, "green", "filled");
scatter(X(C3start:C3end), Y(C3start:C3end), 20, "blue", "filled");
scatter(X(C4start:C4end), Y(C4start:C4end), 20, "magenta", "filled");
% scatter(X(C5start:C5end), Y(C5start:C5end), 20, "black", "filled");
hold off
axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')






