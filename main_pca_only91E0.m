clear all
clc

%% Data import and prep
% INFO: automaticke nastavenie adresara na miesto, kde je tento subor
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

% INFO: az po dalsi komentar INFO je vsetko rovnake ako v subore main_pca.m
% sentinel data
dataSent = readmatrix("112x72.csv");
dataMonoSent = readmatrix("mono_11x72_2018.csv");

% point cloud data
dataRM = readmatrix("RM_curves_rev2.3.csv");
% dataRM = readmatrix("RM_rev2.3_onlyVegHigh.csv");

dataMonoRM = readmatrix("RM_curvesMono_inSVK.csv");
% dataMonoRM = readmatrix("RM_mono_onlyVegHigh.csv");


index = @(i) ((i-1)*4 + 1):1:(i*4);

selected = [];
metrics = [ 1:1:22 ];
for i = metrics
	selected = [ selected index(i) ]; %#ok<AGROW>
end

data = dataSent;
% data = dataRM(:, selected);
% data = [dataSent dataRM(:, selected)];

dataMono = dataMonoSent;
% dataMono = dataMonoRM(:, selected);
% dataMono = [dataMonoSent dataMonoRM(:, selected)];

% INFO: spojenie dat pre biotopy (iba 91E0, prvych 22 riadkov) a dat pre monokultury
dataAll = [data(1:22, :); dataMono];

dataMin = min(dataAll);
dataMax = max(dataAll);
% normovanie
dataScaled = (dataAll - dataMin) ./ (dataMax - dataMin);
dataScaled(isnan(dataScaled)) = 0;

[coef, score, ~, ~, explained, ~] = pca(dataScaled, "Algorithm","eig");
fprintf("explained: %.2f\n", sum(explained(1:3)));

% INFO: znovu len pomocne indexi
n1 = 22;
n2 = 11;
N = n1 + n2;

C1start = 1;
C2start = n1 + 1;

C1end = n1;
C2end = n1 + n2;

ind = [C1start, C1end; C2start, C2end];

X = [score(:,1)];
Y = [score(:,2)];
Z = [score(:,3)];

figure
hold on
scatter(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), 20, "red", "filled");
scatter(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', 'Monokultura', 'Location','southeast')

% 3D graf
% figure("Name","Z = " + num2str(z) + ". súradnica po PCA")
% scatter3(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), Z(ind(1,1):ind(1,2)), 20, "red", "filled");
% hold on
% scatter3(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), Z(ind(2,1):ind(2,2)), 50, "black", "*");
% hold off
% % axis([-0.1 1.2 -0.1 1.2])
% legend('91E0', 'Monokultura', 'Location','southeast')
% title("Z = " + num2str(z) + ". súradnica po PCA")

% INFO: vypocet priemerov a std
means = cell(2,1);
sigmas = cell(2,2);

for i = 1:2
	means{i} = mean( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,1} = means{i} - std( score( ind(i,1):ind(i,2), : ) );
	sigmas{i,2} = means{i} + std( score( ind(i,1):ind(i,2), : ) );
end

figure; 
hold on
features = size(score,2);
% features = ;

plot([0,0],[0,0],'-r');
plot([0,0],[0,0],'-k');

plot(1:1:features, means{1}(1:features),'.-r', "LineWidth", 4, "MarkerSize", 20)
plot(1:1:features, sigmas{1,1}(1:features),':r', "LineWidth", 2)
plot(1:1:features, sigmas{1,2}(1:features),':r', "LineWidth", 2)

plot(1:1:features, means{2}(1:features),'.-k', "LineWidth", 3, "MarkerSize", 20)
plot(1:1:features, sigmas{2,1}(1:features),':k', "LineWidth", 1.5)
plot(1:1:features, sigmas{2,2}(1:features),':k', "LineWidth", 1.5)

legend('91E0', 'Monokultura', 'Location','southeast')
hold off






