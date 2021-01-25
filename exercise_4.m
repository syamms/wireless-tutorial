clear; close all; clc; setup;

nSamples = 1e6;

% * Number of i.i.d. variables
nVariables = 2.^(0 : 5);
nCases = length(nVariables);

% * CSCG distribution [h ~ CN(0, 1)]
cscgSample = sqrt(1 / 2) * (randn(nSamples, nVariables(end)) + 1i * randn(nSamples, nVariables(end)));

figure('name', 'PDF of normal distribution');
tiledlayout(nCases, 1, 'tilespacing', 'compact');

for iVariable = 1 : length(nVariables)
	nexttile;
	histogram(sum(abs(cscgSample(:, 1 : nVariables(iVariable))).^2, 2), 'normalization', 'probability');
	xlim([0, inf]);
	legend(sprintf('$n = %d$', nVariables(iVariable)));
	xlabel('Value');
	ylabel('Probability');
end
