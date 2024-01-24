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

NumOfLIDAR = size(dataLIDAR,2)/4; % number of lidar metrics
NumOfSENTI = size(dataSENTI,2)/4; % number of sentinel metrics

% choose metrics
% metricsLIDAR = 1:NumOfLIDAR;
% metricsLIDAR = [5 13 20];
% metricsLIDAR = 1:17;
metricsLIDAR = [];

metricsSENTI  = 1:NumOfSENTI;
% metricsSENTI = [];

% choose statistics
statistics = [1,1,1,1]; % [mean, std, min, max]

selectedLIDAR = statistics'*ismember(1:NumOfLIDAR,metricsLIDAR);
selectedLIDAR = logical(selectedLIDAR(:));
selectedSENTI = statistics'*ismember(1:NumOfSENTI,metricsSENTI);
selectedSENTI = logical(selectedSENTI(:));

data = [dataSENTI(:, selectedSENTI), dataLIDAR(:, selectedLIDAR)];
dataMono = [dataSENTImono(:, selectedSENTI) dataLIDARmono(:, selectedLIDAR)];
data = [data; dataMono];
data = data(selectSamples,:);

% choose if PCA should be used and what % of explained variance to take
PCA = false;
explainedThreshold = 99.9;

%% Remove the variables with low variance
selectedVariables = true(1,size(data,2));
CVThreshold = 1; % coefficient of variation threshold
StdThreshold = 0.001;

% Using total variance
means = mean(data);
stds = std(data);
CV = (stds ./ means) * 100; % Coefficient of variation for each column

selectedVariables = selectedVariables & (stds > StdThreshold);
selectedVariables = selectedVariables & (CV >= CVThreshold); % columns with variance above the threshold

selectedVariables_id = find(selectedVariables);
removedVariables_id = find(~selectedVariables);

data = data(:, selectedVariables); % Select columns with sufficient variance

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

%% Plot data
figure
imagesc(data);
% imagesc((data-mean(data))./std(data)); % standardized
colorbar;
title('2D Plot of all data');

sigma = cov((data-mean(data))./std(data)); % aby to bolo lepsie vidno na grafe
% sigma = cov(data);
% invSigma = inv(sigma);
figure
imagesc(sigma);
yline(7*2,"LineWidth",3)
yline(18*2,"LineWidth",3)
xline(7*2,"LineWidth",3)
xline(18*2,"LineWidth",3)
colorbar;
title('2D Plot of covariance matrx');

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

	[coef, dataAfterPCA, ~, ~, explained, ~] = pca(dataScaled, "Algorithm","eig");

	sum_explained = 0;
	
	for index = 1:length(explained)
		sum_explained = sum_explained + explained(index);

		if sum_explained > explainedThreshold
			break;
		end
	end

	fprintf("Found index %d\nexplained: %.3f\n", index, sum_explained);

	dataAfterPCA = dataAfterPCA(:, 1:index);

end

%% LDA
if (PCA)
	LDAcls = fitcdiscr(dataAfterPCA, class, 'DiscrimType', 'linear');
else
	LDAcls = fitcdiscr(data, class, 'DiscrimType', 'linear');
end

% predClass = classify(data,data,class); % LDA using classify function

%% Predict Classes
if (PCA)
	predClass = predict(LDAcls, dataAfterPCA);
else
	predClass = predict(LDAcls, data);
end


%% Success rate
figure
CM = confusionmat(class, predClass); % Confusion matrix

format bank

successRateClass = diag(CM)./nk' * 100 % Success rate for each class

successRate = sum(diag(CM))/n * 100; % Total Success Rate:
fprintf('Total Success Rate: %.2f %%\n', successRate);

confusionchart(class,predClass)

%% Vypocet a vykreslenie projekcii
classNames = ["91E0", "91F0","91G0","9110","Mono"];
yLineType = ["r:","g:","b:","m:","k:"];

figure("WindowState","maximized","Name","PCA: Lidar");
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
			dotProducts = sum(dataAfterPCA' .* normal)';
		else
			dotProducts = sum(data' .* normal)';
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
		ylim([0.8, 5.2])

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

figure("WindowState","maximized","Name","Projections 2D plots");
tiledlayout(5, 4);

% ciselna kombinacia pre dvojicu (9110, Mono) = (4, 5)
xAxis_i = 4; xAxis_j = 1;

for i = 1:5
	for j = 1:5
		if (i == j)
			continue
		end

		nexttile
		hold on
		% vykreslenie projekcii
		plot(c91E0{xAxis_i, xAxis_j}(:,1), c91E0{i,j}(:,1), 'r.', "MarkerSize", 15)
		plot(c91F0{xAxis_i, xAxis_j}(:,1), c91F0{i,j}(:,1), 'g.', "MarkerSize", 15)
		plot(c91G0{xAxis_i, xAxis_j}(:,1), c91G0{i,j}(:,1), 'b.', "MarkerSize", 15)
		plot(c9110{xAxis_i, xAxis_j}(:,1), c9110{i,j}(:,1), 'm.', "MarkerSize", 15)
		plot(cmono{xAxis_i, xAxis_j}(:,1), cmono{i,j}(:,1), 'k.', "MarkerSize", 15)
		hold off
		xlabel(classNames(xAxis_i) + " ~ " + classNames(xAxis_j), "FontWeight", "bold")
		ylabel(classNames(i) + " ~ " + classNames(j), "FontWeight", "bold")
	end
end

%% CDA
if (PCA)
	X = dataAfterPCA;
else
	X = data;
end

% Within-Class variance W
Xi = [];
W = zeros(size(X,2));
for i = 1:NumOfClasses
	Xi = X(class == i, :);

	W = W + Xi' * Xi;
end

W = W / (n - 1);

% Sample Variance-Covariance matrix V = B + W
V = (X' * X) / (n - 1);

% Between-Class variance B = V - W
B = V - W;

% F-ratio -> B.u = rho*W.u - generalized eigenvalue problem
[u, lambda] = eig(V\B);
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



