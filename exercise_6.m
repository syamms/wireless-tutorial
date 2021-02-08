clear; close all; clc; setup;
nTxs = 2;
nRxs = 2;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[symbol] = mod_qpsk(bitStream, 2 * bitPower);

numDetBer = zeros(length(snrDb), nChannels);
numAlamoutiBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		fading = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
		noise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(symbol)) + 1i * randn(nRxs, length(symbol)));
		[txDetSymbol] = prec_det(symbol, fading);
		[txAlamoutiSymbol] = prec_alamouti(symbol);
		rxDetSymbol = fading * txDetSymbol + noise;
		rxAlamoutiSymbol = fading * txAlamoutiSymbol + noise;
		[deDetSymbol] = comb_det(rxDetSymbol, fading);
		[deAlamoutiSymbol] = comb_alamouti(rxAlamoutiSymbol, fading);
		[detBit] = demod_qpsk(deDetSymbol);
		[alamoutiBit] = demod_qpsk(deAlamoutiSymbol);
		numDetBer(iSnr, iChannel) = sum(xor(bitStream, detBit)) / nBits;
		numAlamoutiBer(iSnr, iChannel) = sum(xor(bitStream, alamoutiBit)) / nBits;
	end
end
numDetBer = mean(numDetBer, 2);
numAlamoutiBer = mean(numAlamoutiBer, 2);

figure('name', 'BER of QPSK over 2-by-2 MIMO i.i.d. Rayleigh fading channel by DET and Alamouti');
semilogy(snrDb, numDetBer);
hold on;
semilogy(snrDb, numAlamoutiBer);
grid on;
legend('DET', 'Alamouti', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');
