clear
clc
close all

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

data = [dataSENTI(:, selectedSENTI), dataLIDAR(:, selectedLIDAR)];
dataMono = [dataSENTImono(:, selectedSENTI) dataLIDARmono(:, selectedLIDAR)];
data = [data; dataMono];
data = data(selectSamples,:);

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
% imagesc(data);
imagesc((data-mean(data))./std(data)); % standardized
colorbar;
title('2D Plot of all data');

sigma = cov((data-mean(data))./std(data));
% invSigma = inv(sigma);
figure
imagesc(sigma);
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

%% LDA
LDAcls = fitcdiscr(data,class, 'DiscrimType', 'linear');

% predClass = classify(data,data,class); % LDA using classify function

%% Predict Classes
predClass = predict(LDAcls,data);

%% Success rate
figure
CM = confusionmat(class,predClass); % Confusion matrix

format bank

successRateClass = diag(CM)./nk' * 100 % Success rate for each class

successRate = sum(diag(CM))/n * 100; % Total Success Rate:
fprintf('Total Success Rate: %.2f %%\n', successRate);

confusionchart(class,predClass)

%% Vypocet a vykreslenie projekcii
classNames = ["91E0", "91F0","91G0","9110","Mono"];
yLineType = ["r:","g:","b:","m:","k:"];

figure("WindowState","maximized");
tiledlayout(5, 4);

for i = 1:5
	for j = 1:5
		if (i == j)
			continue
		end

		a = LDAcls.Coeffs(i,j).Linear;
		c0 = LDAcls.Coeffs(i,j).Const;

		n = a / norm(a);
		DB = -c0 / norm(a);

		% vypocet projekcii
		dotProducts = sum(data' .* n)';

		coords = [dotProducts, class];

		c91E0 = coords(coords(:,2) == 1, :);
		c91F0 = coords(coords(:,2) == 2, :);
		c91G0 = coords(coords(:,2) == 3, :);
		c9110 = coords(coords(:,2) == 4, :);
		cmono = coords(coords(:,2) == 5, :);

		nexttile
		title(classNames(i) + " ~ " + classNames(j), "coeff(" + num2str(i) + "," + num2str(j) + ")")
		hold on
		% vykreslenie projekcii
		plot(c91E0(:,1), c91E0(:,2), 'r.')
		plot(c91F0(:,1), c91F0(:,2), 'g.')
		plot(c91G0(:,1), c91G0(:,2), 'b.')
		plot(c9110(:,1), c9110(:,2), 'm.')
		plot(cmono(:,1), cmono(:,2), 'k.')
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

	end
end

% legend("91E0", "91F0","91G0","9110","Mono")






