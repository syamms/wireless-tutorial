clear; close all; clc; setup;

nSamples = 1e6;
standardDeviationDb = 8;
maxDisplayValue = 3e1;

% * Normal distribution [x ~ N(0, sigma^2)]
normalSample = standardDeviationDb * randn(nSamples, 1);

% * Log-normal distribution [10log10(x) ~ N(0, sigma^2)]
logNormalSample = db2pow(normalSample);

figure('name', 'PDF of normal and log-normal distribution');
tiledlayout(2, 1, 'tilespacing', 'compact');

nexttile;
histogram(normalSample, 'normalization', 'probability');
legend('Normal');
xlabel('Value');
ylabel('Probability');

nexttile;
histogram(logNormalSample, 'normalization', 'probability');
xlim([0, maxDisplayValue]);
legend('Log-Normal');
xlabel('Value');
ylabel('Probability');
