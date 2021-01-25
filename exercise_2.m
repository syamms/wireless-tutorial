clear; close all; clc; setup;

nSamples = 1e6;

% * CSCG distribution [h ~ CN(0, 1)]
cscgSample = sqrt(1 / 2) * (randn(nSamples, 1) + 1i * randn(nSamples, 1));

figure('name', 'PDF of normal distribution');
tiledlayout(3, 1, 'tilespacing', 'compact');

nexttile;
histogram(abs(cscgSample), 'normalization', 'probability');
xlim([0, inf]);
legend('Magnitude');
xlabel('Value');
ylabel('Probability');

nexttile;
histogram(angle(cscgSample), 'normalization', 'probability');
xlim([-pi, pi]);
legend('Phase');
xlabel('Value');
ylabel('Probability');

nexttile;
histogram(abs(cscgSample).^2, 'normalization', 'probability');
xlim([0, inf]);
legend('Squared magnitude');
xlabel('Value');
ylabel('Probability');
