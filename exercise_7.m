clear; close all; clc; setup;
nTxs = 2;
nRxs = 2;
corTx = [0, 0.9];
corRx = 0;
nChannels = 1e6;

strength = zeros(min(nTxs, nRxs), length(corTx), nChannels);

for iCorTx = 1 : length(corTx)
	for iChannel = 1 : nChannels
		fading = fading_kronecker(nTxs, nRxs, corTx(iCorTx), corRx);
		strength(:, iCorTx, iChannel) = svd(fading);
	end
end

figure('name', 'Channel eigenvalue distribution vs correlation coefficient');
tiledlayout(min(nTxs, nRxs), length(corTx), 'tilespacing', 'compact');
for iCorTx = 1 : length(corTx)
	for iStream = 1 : min(nTxs, nRxs)
		nexttile;
		histogram(strength(iStream, iCorTx, :), 'normalization', 'probability');
		xlim([0, inf]);
		legend(sprintf('Eigenvalue $%d $', iStream));
		xlabel('Value');
		ylabel('Probability');
		title(sprintf('$t = %s$', num2str(corTx(iCorTx))))
	end
end
