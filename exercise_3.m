clear; close all; clc; setup;

nSamples = 1e6;

% * Racian factor [K]
racianFactor = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6];
nCases = length(racianFactor);

% * Phase shift [alpha]
phaseShift = 2 * pi * rand;

% * LoS component [\bar{h}]
losComponent = exp(1i * phaseShift);

% * NLoS component [\tilde{h}]
nlosComponent = sqrt(1 / 2) * (randn(nSamples, 1) + 1i * randn(nSamples, 1));

channel = zeros(nSamples, nCases);
for iCase = 1 : nCases
	channel(:, iCase) = sqrt(racianFactor(iCase) / (1 + racianFactor(iCase))) * losComponent + sqrt(1 / (1 + racianFactor(iCase))) * nlosComponent;
end

figure('name', 'PDF of normal distribution');
tiledlayout(nCases, 1, 'tilespacing', 'compact');

for iCase = 1 : nCases
	nexttile;
	histogram(abs(channel(:, iCase)), 'normalization', 'probability');
	xlim([0, inf]);
	legend(sprintf('$K = %d$', racianFactor(iCase)));
	xlabel('Value');
	ylabel('Probability');
end
