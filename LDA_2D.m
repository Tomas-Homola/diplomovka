clc;
close all;
clear;

% rng(1);  % Set the random seed for reproducibility

%% generate data
% Number of samples per class
numSamples = [100,100,100];

% Mean vectors for each class
mu1 = [1 2];
mu2 = [5 5];
mu3 = [8 2];
mu = [mu1;mu2;mu3];

% Covariance matrix (same for all classes)
sigma = [2 0.5; 0.5 1];
invSigma = inv(sigma);

% Generate synthetic data for each class
data1 = mvnrnd(mu1, sigma, numSamples(1));
data2 = mvnrnd(mu2, sigma, numSamples(2));
data3 = mvnrnd(mu3, sigma, numSamples(3));

% Combine data from all classes
allData = [data1; data2; data3];

% Create labels
labels = [ones(numSamples(1), 1); 2 * ones(numSamples(2), 1); 3 * ones(numSamples(3), 1)];

% Specify colors for each class
col = [[0 0.4470 0.7410];
       [0.8500 0.3250 0.0980];
       [0.9290 0.6940 0.1250]];

%% LDA
% Train Linear Discriminant Analysis (LDA) classifier
LDAcls = fitcdiscr(allData, labels, 'DiscrimType', 'linear');

%% 2D Plot
figure(1);
hold on;
% Point colors for scatter plot
colors = [repmat(col(1,:), numSamples(1), 1);
          repmat(col(2,:), numSamples(2), 1);
          repmat(col(3,:), numSamples(3), 1)];

% Scatter plot of the synthetic data
scatter(data1(:, 1), data1(:, 2), 30, col(1,:),'filled');
scatter(data2(:, 1), data2(:, 2), 30, col(2,:),'filled');
scatter(data3(:, 1), data3(:, 2), 30, col(3,:),'filled');
title('Synthetic Data with 3 Classes');
xlabel('Feature 1');
ylabel('Feature 2');
legend();
hold off;

% Construct meshgrid
[X, Y] = meshgrid(linspace(min(allData(:, 1)), max(allData(:, 1)), 100), ...
                      linspace(min(allData(:, 2)), max(allData(:, 2)), 100));
for k = 1:3
    % Plot contours of PDF
    Z_PDF = mvnpdf([X(:) Y(:)], mu(k,:), sigma);
    Z_PDF = reshape(Z_PDF, size(X));
    figure(1);
    hold on
    contour(X, Y, Z_PDF,20,'EdgeColor',col(k,:), 'DisplayName', ['PDF Contours Class ' num2str(k)]);

    % Compute coefficients for linear discriminant (LD) functions (exact)
    const = - 0.5 * mu(k, :) * invSigma * mu(k, :)';
    lin = mu(k, :) * invSigma;
    % quad = - 0.5 * invSigma;
    % Compute LD functions values
    Z_LD = [X(:) Y(:)] * lin' + const;
    Z_LD = reshape(Z_LD, size(X));
    % Z_QD = dot([X(:) Y(:)],[X(:) Y(:)]*quad,2) + [X(:) Y(:)] * lin' + const;
    % Z_QD = reshape(Z_QD, size(X));
    % Plot contours of LD function
    % contour(X, Y, Z_LD,'EdgeColor',col(k,:), 'DisplayName', ['LD Contours, Class ' num2str(k)]);
    
    % Compute (estimate) decision boundary from sample
    x = linspace(min(allData(:, 1)), max(allData(:, 1)), 100);
    for i = k+1:3
        y = (LDAcls.Coeffs(k,i).Const + LDAcls.Coeffs(k,i).Linear(1) * x)/LDAcls.Coeffs(k,i).Linear(2);
        
        % Compute decision boundary function
        Z_DB = [X(:) Y(:)] * LDAcls.Coeffs(k,i).Linear + LDAcls.Coeffs(k,i).Const;
        Z_DB = reshape(Z_DB, size(X));
        
        % Plot estimated decision boundary
        contour(X, Y, Z_DB, [0 0],'--', 'DisplayName', ['Estimated Decision Boundary ' num2str(k),'-', num2str(i)], 'LineColor', 'k');
    end

    % Plot means
    scatter(mu(k, 1), mu(k, 2), 100, 'k', 'x', 'LineWidth', 2, 'DisplayName', ['Mean Class ' num2str(k)]);
    hold off

    figure(2);
    title('PDF and LD functions');
    hold on
    colors = repmat(reshape(col(k,:), 1, 1, []), [size(X), 1]);
    surf(X, Y, 500*Z_PDF,colors);
    surf(X, Y, Z_LD,colors);
    alpha(0.5);
    view(3)
    hold off;
end

%% Predict Classes
predictedLabels = predict(LDAcls,allData);

figure
confusionchart(labels,predictedLabels)