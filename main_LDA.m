clear
clc
% close all

%% Data import and prep
[filePath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(filePath);

selectedClasses = logical([1,1,1,1,1]); % [91E0,91F0,91G0,9110,Monokultura]
nk = [22,28,31,31,11]; % number of samples in classes
selectSamples = repelem(selectedClasses, nk);
nk = nk(selectedClasses);
NumOfClasses = nnz(selectedClasses);
n = sum(nk); % total number of samples
class = repelem(1:NumOfClasses, nk)';

% INFO: nacitanie dat, najprv Sentinel data, potom point cloud data
% sentinel data
dataSENTI = readmatrix("112x72.csv");
dataSENTImono = readmatrix("mono_11x72_2018.csv");

% lidar data
% dataLIDAR = readmatrix("RM_curves_rev2.3.csv");
dataLIDAR = readmatrix("repMetrics\RM_biotops_allVeg.csv");
dataLIDARmono = readmatrix("repMetrics\RM_monoculture_allVeg.csv");

% lidar metrics names
namesLIDAR = ["Hmax", "Hmean","Hmedian","Hp25","Hp75", "Hp95"...
			  "PPR","DAM_z","BR_bellow_1","BR_1_2","BR_2_3","BR_above_3","BR_3_4","BR_4_5","BR_bellow_5","BR_5_20","BR_above_20" ...
			  "Coeff_var_z", "Hkurt", "Hskew", "Hstd", "Hvar"]; %, "Shannon"
namesLIDAR = repelem(namesLIDAR, 4)';

namesReprMetrics = ["_mean", "_std", "_min", "_max"]';
namesReprMetrics = repmat(namesReprMetrics, size(dataLIDAR,2)/4, 1);

namesLIDAR = namesLIDAR + namesReprMetrics;

NumOfLIDAR = size(dataLIDAR,2)/4; % number of lidar metrics
NumOfSENTI = size(dataSENTI,2)/4; % number of sentinel metrics

% choose metrics
metricsLIDAR = 1:NumOfLIDAR;
% metricsLIDAR = [5 13 20];
% metricsLIDAR = 1:17;
% metricsLIDAR = [];

% metricsSENTI  = 1:NumOfSENTI;
metricsSENTI = [];

% choose statistics
statistics = [1,1,1,1]; % [mean, std, min, max]

selectedLIDAR = statistics'*ismember(1:NumOfLIDAR,metricsLIDAR);
selectedLIDAR = logical(selectedLIDAR(:));
selectedSENTI = statistics'*ismember(1:NumOfSENTI,metricsSENTI);
selectedSENTI = logical(selectedSENTI(:));

selectedNamesLIDAR = namesLIDAR(selectedLIDAR);

data = [dataSENTI(:, selectedSENTI), dataLIDAR(:, selectedLIDAR)];
dataMono = [dataSENTImono(:, selectedSENTI) dataLIDARmono(:, selectedLIDAR)];
data = [data; dataMono];
data = data(selectSamples,:);

% choose if PCA should be used and what % of explained variance to take
PCA = true;
explainedThreshold = 99.9;

%% ANOVA
pValues = zeros(1, size(data,2));
for i = 1:size(data,2)
	[pValues(i), ~] = anova1(data(:, i), class, "off");
end

% A = anova(class, data(:, nn));
% [pValue, ~] = anova1(data(:, nn), class, "off");
x = 1:1:size(data,2);

alpha = 0.05;
% The small p-value indicates that differences between column means are significant.
passedANOVA = pValues < alpha; 
figure
plot(x(pValues < alpha), pValues(pValues < alpha), '.g', 'MarkerSize', 15)
hold on
plot(x(pValues >= alpha), pValues(pValues >= alpha), '.r', 'MarkerSize', 15)
for i = 1:length(x) 
	plot([x(i) x(i)], [0 pValues(i)], ':k')
end
xline(6*nnz(statistics)+0.5, "--k","LineWidth",1)
xline(17*nnz(statistics)+0.5,"--k","LineWidth",1)
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')

%% KS-test
pValuesKS = zeros(1, size(data,2));
for i = 1:size(data,2)
	[~, pValuesKS(i)] = kstest((data(:, i) - mean(data(:, i))) ./ std(data(:, i)));
end

x = 1:1:size(data,2);

alpha = 0.05;
passedKS = pValuesKS > alpha; % Small values of p cast doubt on the validity of the null hypothesis.
figure
plot(x(pValuesKS < alpha), pValuesKS(pValuesKS < alpha), '.r', 'MarkerSize', 15)
hold on
plot(x(pValuesKS >= alpha), pValuesKS(pValuesKS >= alpha), '.g', 'MarkerSize', 15)
for i = 1:length(x) 
	plot([x(i) x(i)], [0 pValuesKS(i)], ':k')
end
xline(6*nnz(statistics)+0.5, "--k","LineWidth",1)
xline(17*nnz(statistics)+0.5,"--k","LineWidth",1)
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')

%% Remove the variables with low variance
selectedVariables = true(1,size(data,2));
CVThreshold = 1; % coefficient of variation threshold
StdThreshold = 0.001;

% Using total variance
means = mean(data);
stds = std(data);
CV = (stds ./ means) * 100; % Coefficient of variation for each column

% selectedVariables = selectedVariables & (stds > StdThreshold);
% selectedVariables = selectedVariables & (CV >= CVThreshold); % columns with variance above the threshold
selectedVariables = selectedVariables & ~passedKS;

selectedVariables_id = find(selectedVariables);
removedVariables_id = find(~selectedVariables);

fprintf("Removed variables:\n");
selectedNamesLIDAR(~selectedVariables)

data = data(:, selectedVariables); % Select columns with sufficient variance
selectedNamesLIDAR = selectedNamesLIDAR(selectedVariables);

% Using within class variance 
% selectedVariables = true(1,nnz(selectedVariables));
% for k = 1:NumOfClasses
%     classData = data(class == k, :);
%     classMean = mean(classData);
%     classStd = std(classData);
%     withinClassCV = (classMean ./ classStd) * 100;
%     selectedVariables = selectedVariables & (classStd > StdThreshold); % unselect columns with zero std
%     selectedVariables = selectedVariables & (withinClassCV > CVThreshold); % unselect column with low CV
% end
% selectedVariables_id = find(selectedVariables);
% removed_id = find(~selectedVariables);
% 
% data = data(:,selectedVariables);

%% Plot data correlation
% figure
% imagesc(data);
% % imagesc((data-mean(data))./std(data)); % standardized
% colorbar;
% title('2D Plot of all data');

sigma = cov((data-mean(data))./std(data)); % aby to bolo lepsie vidno na grafe
% sigma = corr(data);
% sigma = cov(data);
% invSigma = inv(sigma);
figure
imagesc(sigma);
yline(6*nnz(statistics)+0.5, "-k","LineWidth",3)
yline(17*nnz(statistics)+0.5,"-k","LineWidth",3)
xline(6*nnz(statistics)+0.5, "-k","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-k","LineWidth",3)
colorbar;
title('2D Plot of covariance matrx');
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
yticks(1:1:size(data,2))
yticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
colormap("turbo")

% Plot within class covariance matrices
% for k = 1:NumOfClasses
%     classData = data(class == k, :);
%     figure
%     imagesc(classData);
%     colorbar;
%     title(['2D Plot of classData ',num2str(k)]);
% 
%     sigma = cov(classData);
%     invSigma = inv(sigma);
%     figure
%     imagesc(sigma);
%     colorbar;
%     title(['2D Plot of covariance matrx ',num2str(k)]);
% end

%% PCA
if (PCA)
	dataScaled = [];
	dataScaled = (data - mean(data)) ./ std(data);
	% 	dataScaled = zscore(data); % funkcia na standardizaciu

	[coef, dataAfterPCA, variances, ~, explained, ~] = pca(dataScaled, "Algorithm","eig");

	index = find(cumsum(explained) > explainedThreshold, 1);

	fprintf("Found index %d\nexplained: %.3f\n", index, sum(explained(1:index)));

	dataAfterPCA = dataAfterPCA(:, 1:index);

end

%% Biplot pre PCA
statisticsPCA = [1,1,1,1]; % [mean, std, min, max]

plotCoef = statisticsPCA'*ismember(1:NumOfLIDAR,metricsLIDAR);
plotCoef = logical(plotCoef(:));
plotCoef = plotCoef(selectedVariables);

figure
% biplot(coef(:,1:2),"VarLabels",selectedNamesLIDAR)
biplotG(coef(plotCoef, :), [], "VarLabels",selectedNamesLIDAR(plotCoef))

maxlen = max(sqrt(coef(:,1).^2 + coef(:,2).^2));
scores = dataAfterPCA./max(max(abs(dataAfterPCA)))*maxlen;
hold on
scatter3(scores(class == 1, 1), scores(class == 1, 2), scores(class == 1, 3),"or",'MarkerFaceColor', 'r')
scatter3(scores(class == 2, 1), scores(class == 2, 2), scores(class == 2, 3),"og",'MarkerFaceColor', 'g')
scatter3(scores(class == 3, 1), scores(class == 3, 2), scores(class == 3, 3),"ob",'MarkerFaceColor', 'b')
scatter3(scores(class == 4, 1), scores(class == 4, 2), scores(class == 4, 3),"om",'MarkerFaceColor', 'm')
scatter3(scores(class == 5, 1), scores(class == 5, 2), scores(class == 5, 3),"ok",'MarkerFaceColor', 'k')
hold off
% view(3)

%% Prispevky prediktorov po PCA #1
coef2 = normalize(abs(coef), 'norm',1);

% figure
% biplot(coef(:,1:2),"VarLabels",selectedNamesLIDAR);

figure
plotBarPCA(coef, 1)
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
xline(6*nnz(statistics)+0.5, "-g","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-g","LineWidth",3)

figure
plotBarPCA(coef, 2)
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
xline(6*nnz(statistics)+0.5, "-g","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-g","LineWidth",3)

figure
plotBarPCA(coef, 3)
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
xline(6*nnz(statistics)+0.5, "-g","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-g","LineWidth",3)

%% Cross-validation
rep = 200;
accuracies = zeros(rep, 1);
k = 5;

for i = 1:rep
	cv = cvpartition(n, "KFold", k);
	if (PCA)
		LDAmodel = fitcdiscr(dataAfterPCA, class, 'CVPartition', cv, 'DiscrimType', 'linear');
	else
		LDAmodel = fitcdiscr(data, class, 'CVPartition', cv, 'DiscrimType', 'linear');
	end
	error = kfoldLoss(LDAmodel);
	accuracies(i) = (1 - error) * 100;
end

meanAccuracy = mean(accuracies);
stdAccuracy = std(accuracies);

fprintf("CV result accuracy: %.2f %% +- %.2f %%\n", meanAccuracy, stdAccuracy);

%% LDA - fit model
if (PCA)
	LDAcls = fitcdiscr(dataAfterPCA, class, 'DiscrimType', 'linear');
else
	LDAcls = fitcdiscr(data, class, 'DiscrimType', 'linear');
end

% predClass = classify(data,data,class); % LDA using classify function

%% LDA - Predict Classes
if (PCA)
	predClass = predict(LDAcls, dataAfterPCA);
else
	predClass = predict(LDAcls, data);
end


%% LDA - Success rate
figure
CM = confusionmat(class, predClass); % Confusion matrix

format bank

successRateClass = diag(CM)./nk' * 100 % Success rate for each class

successRate = sum(diag(CM))/n * 100; % Total Success Rate:
fprintf('Total Success Rate: %.2f %%\n', successRate);

confusionchart(class,predClass)

%% LDA - Vypocet a vykreslenie projekcii
classNames = ["91E0", "91F0","91G0","9110","Mono"];
yLineType = ["r:","g:","b:","m:","k:"];

figure("WindowState","maximized","Name","PCA: Lidar+Sentinel");
tiledlayout(5, 4);

c91E0 = cell(5, 4);
c91F0 = cell(5, 4);
c91G0 = cell(5, 4);
c9110 = cell(5, 4);
cmono = cell(5, 4);

for i = 1:5
	for j = 1:5
		if (i == j)
			continue
		end

		a = LDAcls.Coeffs(i,j).Linear;
		c0 = LDAcls.Coeffs(i,j).Const;

		normal = a / norm(a);
		DB = -c0 / norm(a);

		% vypocet projekcii
		if (PCA)
			dotProducts = dataAfterPCA * normal;
		else
			dotProducts = data * normal;
		end

		coords = [dotProducts, class];

		c91E0{i,j} = coords(coords(:,2) == 1, :);
		c91F0{i,j} = coords(coords(:,2) == 2, :);
		c91G0{i,j} = coords(coords(:,2) == 3, :);
		c9110{i,j} = coords(coords(:,2) == 4, :);
		cmono{i,j} = coords(coords(:,2) == 5, :);

		nexttile
		title(classNames(i) + " ~ " + classNames(j), "coeff(" + num2str(i) + "," + num2str(j) + ")")
		hold on
		% vykreslenie projekcii
		plot(c91E0{i,j}(:,1), c91E0{i,j}(:,2), 'r.', "MarkerSize", 10)
		plot(c91F0{i,j}(:,1), c91F0{i,j}(:,2), 'g.', "MarkerSize", 10)
		plot(c91G0{i,j}(:,1), c91G0{i,j}(:,2), 'b.', "MarkerSize", 10)
		plot(c9110{i,j}(:,1), c9110{i,j}(:,2), 'm.', "MarkerSize", 10)
		plot(cmono{i,j}(:,1), cmono{i,j}(:,2), 'k.', "MarkerSize", 10)
		% horizontalne ciary na zvyraznenie tried
		yline(i, yLineType(i))
		yline(j, yLineType(j))
		% decision boundary
		xl = xline(DB,'--','DB');
		xl.LabelVerticalAlignment = "middle";
		xl.LabelHorizontalAlignment = "center";
		hold off
		% premenovanie tickov na y osi kvoli prehladnosti
		yticks([1 2 3 4 5])
		yticklabels({'91E0', '91F0','91G0','9110','Mono'})
		ylim([0.6, 5.4])

	end
end

% legend("91E0", "91F0","91G0","9110","Mono")

%% Vykreslenie najlepsich projekcii
% ked sa spravi PCA:
% -> iba lidar data: skoro vsetky oddelene, az na 91G0~9110
% -> iba sent. data: vsetko oddelene, najlepsie 9110~Mono a 91G0~Mono
% -> lidar+sentinel: vsetko este lepsie oddelene, najlepsie je dvojica (9110, Mono)

% ked sa NEspravi PCA:
% -> iba lidar data: vsetko oddelene, najlepsie 9110~91E0
% -> iba sent. data: vsetko oddelene, najlepsie 9110~Mono
% -> lidar+sentinel: zhoda, ze to nefunguje dobre kvoli malemu datasetu

% figure("WindowState","maximized","Name","Projections 2D plots");
% tiledlayout(5, 4);
% 
% % ciselna kombinacia pre dvojicu (9110, Mono) = (4, 5)
% xAxis_i = 4; xAxis_j = 5;
% 
% for i = 1:5
% 	for j = 1:5
% 		if (i == j)
% 			continue
% 		end
% 
% 		nexttile
% 		hold on
% 		% vykreslenie projekcii
% 		plot(c91E0{xAxis_i, xAxis_j}(:,1), c91E0{i,j}(:,1), 'r.', "MarkerSize", 15)
% 		plot(c91F0{xAxis_i, xAxis_j}(:,1), c91F0{i,j}(:,1), 'g.', "MarkerSize", 15)
% 		plot(c91G0{xAxis_i, xAxis_j}(:,1), c91G0{i,j}(:,1), 'b.', "MarkerSize", 15)
% 		plot(c9110{xAxis_i, xAxis_j}(:,1), c9110{i,j}(:,1), 'm.', "MarkerSize", 15)
% 		plot(cmono{xAxis_i, xAxis_j}(:,1), cmono{i,j}(:,1), 'k.', "MarkerSize", 15)
% 		hold off
% 		xlabel(classNames(xAxis_i) + " ~ " + classNames(xAxis_j), "FontWeight", "bold")
% 		ylabel(classNames(i) + " ~ " + classNames(j), "FontWeight", "bold")
% 	end
% end

%% CDA - vypocet matic V, W a B + vykreslenie
if (PCA)
	X = dataAfterPCA;
else
	X = data;
end
% X = data(:, passedANOVA);
% nnz(passedANOVA)
% Within-Class variance W
Xi = [];
W = zeros(size(X,2));
for i = 1:NumOfClasses
	Xi = X(class == i, :);

	Xi = Xi - mean(Xi); % treba mat centrovane

	W = W + Xi' * Xi;
end

W = W / (n - 1);

% Sample Variance-Covariance matrix V = B + W
V = cov(X); % aj tu centrovane

% Xc = X - mean(X);
% V2 = (Xc' * Xc) / (n - 1);

% Between-Class variance B = V - W
B = V - W;

% F-ratio -> B.u = rho*W.u - generalized eigenvalue problem
[u, lambda] = eig(B,W);
u = u ./ vecnorm(u);
norms = vecnorm(u);

% zoradenie vlastnych hodnot zostupne + usporiadanie vlastnych vektorov
% podla vlastnych cisel
[lambda_sorted, ind] = sort(diag(lambda),"descend");
u_sorted = u(:, ind);

Xnew = X * u_sorted;

X = Xnew(:,1);
Y = Xnew(:,2);
Z = Xnew(:,3);

figure
hold on
scatter(X(1:22),  Y(1:22), 20, "red", "filled");
scatter(X(23:50), Y(23:50), 20, "green", "filled");
scatter(X(51:81), Y(51:81), 20, "blue", "filled");
scatter(X(82:112), Y(82:112), 20, "magenta", "filled");
scatter(X(113:123), Y(113:123), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')

figure
title 3D
hold on
scatter3(X(1:22),  Y(1:22), Z(1:22), 20, "red", "filled");
scatter3(X(23:50), Y(23:50), Z(23:50), 20, "green", "filled");
scatter3(X(51:81), Y(51:81), Z(51:81), 20, "blue", "filled");
scatter3(X(82:112), Y(82:112),Z(82:112),  20, "magenta", "filled");
scatter3(X(113:123), Y(113:123),Z(113:123), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
view(3)
legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')

%% Prispevky jednotlivych prediktorov #2
% ževraj na znamienku prispevku nezáleží, stačí sa pozerať len na jeho veľkosť
prispevky = coef(:, 1:index) * u_sorted;
prispevky_normed = normalize(abs(prispevky), 'norm', 1);

figure
biplot(prispevky(:, 1:2), "VarLabels", selectedNamesLIDAR);

figure
nn = 1;
bar(prispevky_normed(:,nn) * 100);
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
xline(6*nnz(statistics)+0.5, "-r","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-r","LineWidth",3)
title("CDA " + num2str(nn))

figure
nn = 2;
bar(prispevky_normed(:,nn) * 100);
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
xline(6*nnz(statistics)+0.5, "-r","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-r","LineWidth",3)
title("CDA " + num2str(nn))

%% Biplot pre CDA
statistics = [1,0,0,0]; % [mean, std, min, max]

selectedLIDARpr = statistics'*ismember(1:NumOfLIDAR,metricsLIDAR);
selectedLIDARpr = logical(selectedLIDARpr(:));
selectedLIDARpr = selectedLIDARpr(selectedVariables);

Format = { {'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'};...
		   {'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'};...
		   {'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'};...
		   {'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm'};...
		   {'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'}};

biplotG(prispevky(selectedLIDARpr,:), Xnew,...
	'Groups', class,...
	'VarLabels', selectedNamesLIDAR(selectedLIDARpr),...
	'Format', Format)

%% CDA classification
dsq = zeros(n,NumOfClasses);
invW = inv(W);
for k = 1:NumOfClasses
    dataCenteredClass = data - classMean(k,:); % TODO: dorobit centroidy
    dsq(:,k) = dot(dataCenteredClass,(invW*dataCenteredClass')',2);
end

[~,predClass] = min(dsq,[],2);

% Success rate
figure
CM = confusionmat(class,predClass); % Confusion matrix

format bank

successRateClass = diag(CM)./nk' * 100 % Success rate for each class

successRate = sum(diag(CM))/n * 100; % Total Success Rate:
fprintf('Total Success Rate: %.2f %%\n', successRate);


%%
figure
% imagesc(u_sorted(:,1:4))
imagesc((u_sorted(:,1:4) - mean(u_sorted(:,1:4)) ./ std(u_sorted(:,1:4))))
colorbar
hold on
for i = 1:21
	yline(i*4 + 1,"LineWidth",1)
end
hold off
% yline(7*4,"LineWidth",3)
% yline(18*4,"LineWidth",3)
% xline(7*2,"LineWidth",3)
% xline(18*2,"LineWidth",3)

%% Zistenie uhlov pre normaly z LDA
dots = zeros(5,5);% 3,4
n1 = LDAcls.Coeffs(3,4).Linear;
n1 = n1 / norm(n1);
for i = 1:5
	for j = 1:5
		if (i == j)
			continue
		end
		n2 = LDAcls.Coeffs(i,j).Linear;
		n2 = n2 / norm(n2);
		dots(i,j) = dot(n1, n2);
	end
end

dots

%% LDA - vykreslenie cez projekcie na 2 vybrane normaly 
% 2,5 ~ 2,4, 3,4 ~ 2,4
figure("WindowState","maximized","Name","Projections 2D plots");
tiledlayout(5, 4);
ind1 = [3 4];

for i = 1:5
	for j = 1:5
		if (i == j)
			continue
		end

		newCoords = zeros(n, 2);
		ind2 = [i j];
		for k = 1:n
			newCoords(k,:) = proj(dataAfterPCA(k,:), ind1, ind2, LDAcls);
		end

		X = newCoords(:,1);
		Y = newCoords(:,2);

		nexttile
		hold on
		scatter(X(1:22),  Y(1:22), 20, "red", "filled");
		scatter(X(23:50), Y(23:50), 20, "green", "filled");
		scatter(X(51:81), Y(51:81), 20, "blue", "filled");
		scatter(X(82:112), Y(82:112), 20, "magenta", "filled");
		scatter(X(113:123), Y(113:123), 50, "black", "*");
		hold off
		xlabel('3~4')
		ylabel(num2str(i) + "~" + num2str(j))
		% axis([-0.1 1.2 -0.1 1.2])
% 		legend('91E0', '91F0', '91G0', '9110', 'Monokultura', 'Location','southeast')
	end
end

%% Vykreslenie priemerov pre prediktory
temp = zeros(5, size(data,2));

for i = 1:5
	temp(i, :) = mean(data(class == i, :));
end

temp2 = abs(temp) ./ max(abs(temp));
% temp2 = normalize(temp, 'range');

figure
imagesc(temp2);
xline(6*nnz(statistics)+0.5, "-k","LineWidth",3)
xline(17*nnz(statistics)+0.5,"-k","LineWidth",3)
colorbar;
title('2D Plot of covariance matrx');
xticks(1:1:size(data,2))
xticklabels(selectedNamesLIDAR)
set(gca,'TickLabelInterpreter','none')
colormap("turbo")

%% Pomocne funkcie
function out = proj(x, ind1, ind2, model)
	arguments
		x
		ind1
		ind2
		model ClassificationDiscriminant
	end
	
	out = zeros(2,1);
	n1 = model.Coeffs(ind1(1), ind1(2)).Linear;
	n2 = model.Coeffs(ind2(1), ind2(2)).Linear;

	n1 = n1 / norm(n1);
	n2 = n2 / norm(n2);

	X1 = dot(x, n1); X2 = dot(x, n2);
	N12 = dot(n1, n2);

% 	fprintf("n1.n2: %.4f\n", dot(n1, n2));

	out(1) = ( X1 - X2 * N12 );
	out(2) = ( X2 - X1 * N12 );

	out = out' / (1 - N12 * N12);

end

function plotBarPCA(coef, PC)
	coef2 = normalize(abs(coef), 'norm',1);

	% kod k farbam pre bar plot
	signs = sign(coef(:, PC));
	signs(signs == -1) = 0;
	colors = zeros(size(coef2,1), 3);
	for i = 1:size(coef2,1)
		if signs(i)
			colors(i,:) = [1.0 0.0 0.0];
		else
			colors(i,:) = [0.0 0.0 1.0];
		end
	end

	bar(coef2(:, PC) * 100, "FaceColor","flat","CData", colors);
	title("PC " + num2str(PC))
end




