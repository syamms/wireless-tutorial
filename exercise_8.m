clear; close all; clc; setup;
nTxs = 2;
nRxs = 2;
bitPower = 1;
snrDb = 0 : 5 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
corTx = 0 : 0.2 : 1;
corRx = 0 : 0.2 : 1;
bitStream = round(rand(1, nBits));

[symbol] = mod_qpsk(bitStream, 2 * bitPower);

numDetBer = zeros(length(corTx), length(corRx), length(snrDb), nChannels);
numAlamoutiBer = zeros(length(corTx), length(corRx), length(snrDb), nChannels);

for iCorTx = 1 : length(corTx)
	for iCorRx = 1 : length(corRx)
		for iSnr = 1 : length(snrDb)
			for iChannel = 1 : nChannels
				fading = fading_kronecker(nTxs, nRxs, corTx(iCorTx), corRx(iCorRx));
				noise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(symbol)) + 1i * randn(nRxs, length(symbol)));
				[txDetSymbol] = prec_det(symbol, fading);
				[txAlamoutiSymbol] = prec_alamouti(symbol);
				rxDetSymbol = fading * txDetSymbol + noise;
				rxAlamoutiSymbol = fading * txAlamoutiSymbol + noise;
				[deDetSymbol] = comb_det(rxDetSymbol, fading);
				[deAlamoutiSymbol] = comb_alamouti(rxAlamoutiSymbol, fading);
				[detBit] = demod_qpsk(deDetSymbol);
				[alamoutiBit] = demod_qpsk(deAlamoutiSymbol);
				numDetBer(iCorTx, iCorRx, iSnr, iChannel) = sum(xor(bitStream, detBit)) / nBits;
				numAlamoutiBer(iCorTx, iCorRx, iSnr, iChannel) = sum(xor(bitStream, alamoutiBit)) / nBits;
			end
		end
	end
end
numDetBer = mean(numDetBer, 4);
numAlamoutiBer = mean(numAlamoutiBer, 4);

for iSnr = 1 : length(snrDb)
	figure('name', 'BER of QPSK over 2-by-2 MIMO spatially correlated Rayleigh fading channel by DET and Alamouti');
	surf(corTx, corRx, numDetBer(:, :, iSnr), 'facecolor', 'b');
	hold on;
	surf(corTx, corRx, numAlamoutiBer(:, :, iSnr), 'facecolor', 'r');
	grid on;
	legend('DET', 'Alamouti', 'location', 'sw');
	xlabel('Transmit correlation coefficient');
	ylabel('Receive correlation coefficient');
	zlabel('BER');
	title(sprintf('$SNR = %s$ dB', num2str(snrDb(iSnr))));
end
