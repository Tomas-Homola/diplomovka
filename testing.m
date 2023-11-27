%% Data import and prep
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

% sentinel data
dataSent = readmatrix("112x72.csv");
dataMonoSent = readmatrix("sent_mono_11x72.csv");

% point cloud data
% dataRM = readmatrix("RM_curves_rev2.3.csv");
dataRM = readmatrix("RM_rev2.3_onlyVegHigh.csv");
% dataMonoRM = readmatrix("RM_curvesMono_inSVK.csv");
dataMonoRM = readmatrix("RM_mono_onlyVegHigh.csv");

% zaujimave kombinacie (iba vysoka vegetacia)
% 1) 7 16 18 (78.91 %)
% 2) 5 13 20 -> v 3D (91.2 %)
% 3) 7 18 20 -> v 3D (83.53 %)
% 4) 2 5 13 20 -> v 3D (92.12 %)
index = @(i) ((i-1)*4 + 1):1:(i*4);

selected = [];
metrics = [5 8 13 20];
for i = metrics
	selected = [ selected index(i) ]; %#ok<AGROW> 
end

% data = dataSent;
data = dataRM(:, selected);
% data = [dataSent dataRM(:, selected)];

% dataMono = dataMonoSent;
dataMono = dataMonoRM(:, selected);
% dataMono = [dataMonoSent dataMonoRM(:, selected)];

dataScaled = (data - mean(data))./std(data);
dataScaled(isnan(dataScaled)) = 0;
dataScaled = normalize(dataScaled,"range");

dataMonoScaled = (dataMono - mean(dataMono))./std(dataMono);
dataMonoScaled(isnan(dataMonoScaled)) = 0;
dataMonoScaled = normalize(dataMonoScaled,"range");

[coef, score, ~, ~, explained, ~] = pca(dataScaled, "Algorithm","eig");
fprintf("explained: %.2f\n", sum(explained(1:3)));

scoreMono = dataMonoScaled * coef;

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

X = [score(:,1); scoreMono(:,1)];
Y = [score(:,2); scoreMono(:,2)];
Z = [score(:,3); scoreMono(:,3)];

% X = score(:,1);
% Y = score(:,2);

figure
scatter(X(C1start:C1end), Y(C1start:C1end), 15, "red", "filled");
hold on
scatter(X(C2start:C2end), Y(C2start:C2end), 15, "green", "filled");
scatter(X(C3start:C3end), Y(C3start:C3end), 15, "blue", "filled");
scatter(X(C4start:C4end), Y(C4start:C4end), 15, "magenta", "filled");
scatter(X(C5start:C5end), Y(C5start:C5end), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
title(num2str(metrics))


figure
scatter3(X(C1start:C1end), Y(C1start:C1end), Z(C1start:C1end), 15, "red", "filled");
hold on
scatter3(X(C2start:C2end), Y(C2start:C2end), Z(C2start:C2end), 15, "green", "filled");
scatter3(X(C3start:C3end), Y(C3start:C3end), Z(C3start:C3end), 15, "blue", "filled");
scatter3(X(C4start:C4end), Y(C4start:C4end), Z(C4start:C4end), 15, "magenta", "filled");
scatter3(X(C5start:C5end), Y(C5start:C5end), Z(C5start:C5end), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
title(num2str(metrics))


