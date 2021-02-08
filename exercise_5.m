%% * BPSK/QPSK over AWGN channel
clear; close all; clc; setup;
nTxs = 1;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

anaBer = qfunc(sqrt(2 * db2pow(snrDb)));
numBpskBer = zeros(length(snrDb), nChannels);
numQpskBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		bpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		qpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
		rxBpskSymbol = bpskSymbol + bpskNoise;
		rxQpskSymbol = qpskSymbol + qpskNoise;
		[bpskBit] = demod_bpsk(rxBpskSymbol);
		[qpskBit] = demod_qpsk(rxQpskSymbol);
		numBpskBer(iSnr, iChannel) = sum(xor(bitStream, bpskBit)) / nBits;
		numQpskBer(iSnr, iChannel) = sum(xor(bitStream, qpskBit)) / nBits;
	end
end
numBpskBer = mean(numBpskBer, 2);
numQpskBer = mean(numQpskBer, 2);

figure('name', 'BER of BPSK/QPSK over an AWGN channel');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');

numAwgnBpskBer = numBpskBer;
numAwgnQpskBer = numQpskBer;
clearvars -except numAwgn*

%% * BPSK/QPSK over SISO Rayleigh fading channel
nTxs = 1;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

anaBer = 1 / 2 * (1 - sqrt(db2pow(snrDb) ./ (1 + db2pow(snrDb))));
anaBerUb = 1 ./ (4 * db2pow(snrDb));
numBpskBer = zeros(length(snrDb), nChannels);
numQpskBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		fading = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
		bpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		qpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
		rxBpskSymbol = fading * bpskSymbol + bpskNoise;
		rxQpskSymbol = fading * qpskSymbol + qpskNoise;
		deBpskSymbol = comb_zf(rxBpskSymbol, fading);
		deQpskSymbol = comb_zf(rxQpskSymbol, fading);
		[bpskBit] = demod_bpsk(deBpskSymbol);
		[qpskBit] = demod_qpsk(deQpskSymbol);
		numBpskBer(iSnr, iChannel) = sum(xor(bitStream, bpskBit)) / nBits;
		numQpskBer(iSnr, iChannel) = sum(xor(bitStream, qpskBit)) / nBits;
	end
end
numBpskBer = mean(numBpskBer, 2);
numQpskBer = mean(numQpskBer, 2);

figure('name', 'BER of BPSK/QPSK over SISO Rayleigh fading channel');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, anaBerUb);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');

numZfBpskBer = numBpskBer;
numZfQpskBer = numQpskBer;
clearvars -except numAwgn* numZf*

%% * BPSK/QPSK over SIMO i.i.d. Rayleigh fading channel with 2-rx MRC combining
nTxs = 1;
nRxs = 2;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

gamma = 1 / 2 - 1 / 2 * (1 + 1 ./ db2pow(snrDb)) .^ (- 1 / 2);
anaBer = gamma .^ 2 .* (1 + 2 * (1 - gamma));
anaBerUb = (4 * db2pow(snrDb)) .^ (- nRxs) .* nchoosek(2 * nRxs - 1, nRxs);
numBpskBer = zeros(length(snrDb), nChannels);
numQpskBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		fading = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
		bpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		qpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
		rxBpskSymbol = fading * bpskSymbol + bpskNoise;
		rxQpskSymbol = fading * qpskSymbol + qpskNoise;
		deBpskSymbol = comb_mrc(rxBpskSymbol, fading);
		deQpskSymbol = comb_mrc(rxQpskSymbol, fading);
		[bpskBit] = demod_bpsk(deBpskSymbol);
		[qpskBit] = demod_qpsk(deQpskSymbol);
		numBpskBer(iSnr, iChannel) = sum(xor(bitStream, bpskBit)) / nBits;
		numQpskBer(iSnr, iChannel) = sum(xor(bitStream, qpskBit)) / nBits;
	end
end
numBpskBer = mean(numBpskBer, 2);
numQpskBer = mean(numQpskBer, 2);

figure('name', 'BER of BPSK/QPSK over SIMO i.i.d. Rayleigh fading channel with 2-rx MRC combining');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, anaBerUb);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');

numMrcBpskBer = numBpskBer;
numMrcQpskBer = numQpskBer;
clearvars -except numAwgn* numZf* numMrc*

%% * BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx MRT beamforming
nTxs = 2;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

gamma = 1 / 2 - 1 / 2 * (1 + 1 ./ db2pow(snrDb)) .^ (- 1 / 2);
anaBer = gamma .^ 2 .* (1 + 2 * (1 - gamma));
anaBerUb = (4 * db2pow(snrDb)) .^ (- nTxs) .* nchoosek(2 * nTxs - 1, nTxs);
numBpskBer = zeros(length(snrDb), nChannels);
numQpskBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		fading = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
		bpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		qpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
		[txBpskSymbol] = prec_mrt(bpskSymbol, fading);
		[txQpskSymbol] = prec_mrt(qpskSymbol, fading);
		rxBpskSymbol = fading * txBpskSymbol + bpskNoise;
		rxQpskSymbol = fading * txQpskSymbol + qpskNoise;
		[bpskBit] = demod_bpsk(rxBpskSymbol);
		[qpskBit] = demod_qpsk(rxQpskSymbol);
		numBpskBer(iSnr, iChannel) = sum(xor(bitStream, bpskBit)) / nBits;
		numQpskBer(iSnr, iChannel) = sum(xor(bitStream, qpskBit)) / nBits;
	end
end
numBpskBer = mean(numBpskBer, 2);
numQpskBer = mean(numQpskBer, 2);

figure('name', 'BER of BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx MRT beamforming');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, anaBerUb);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');

numMrtBpskBer = numBpskBer;
numMrtQpskBer = numQpskBer;
clearvars -except numAwgn* numZf* numMrc* numMrt*

%% * Alamouti scheme with BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx
nTxs = 2;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e5;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

gamma = 1 / 2 - 1 / 2 * (1 + 2 ./ db2pow(snrDb)) .^ (- 1 / 2);
anaBer = gamma .^ 2 .* (1 + 2 * (1 - gamma));
anaBerUb = (db2pow(snrDb) * (2 * sqrt(bitPower)) ^ 2 / 8) .^ (- 2);
numBpskBer = zeros(length(snrDb), nChannels);
numQpskBer = zeros(length(snrDb), nChannels);

for iSnr = 1 : length(snrDb)
	for iChannel = 1 : nChannels
		fading = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
		bpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		qpskNoise = sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
		[txBpskSymbol] = prec_alamouti(bpskSymbol);
		[txQpskSymbol] = prec_alamouti(qpskSymbol);
		rxBpskSymbol = fading * txBpskSymbol + bpskNoise;
		rxQpskSymbol = fading * txQpskSymbol + qpskNoise;
		[deBpskSymbol] = comb_alamouti(rxBpskSymbol, fading);
		[deQpskSymbol] = comb_alamouti(rxQpskSymbol, fading);
		[bpskBit] = demod_bpsk(deBpskSymbol);
		[qpskBit] = demod_qpsk(deQpskSymbol);
		numBpskBer(iSnr, iChannel) = sum(xor(bitStream, bpskBit)) / nBits;
		numQpskBer(iSnr, iChannel) = sum(xor(bitStream, qpskBit)) / nBits;
	end
end
numBpskBer = mean(numBpskBer, 2);
numQpskBer = mean(numQpskBer, 2);

figure('name', 'BER of BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, anaBerUb);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');

numAlamoutiBpskBer = numBpskBer;
numAlamoutiQpskBer = numQpskBer;
clearvars -except numAwgn* numZf* numMrc* numMrt* numAlamouti* snrDb

%% * BER comparison
figure('name', 'BER comparison');
semilogy(snrDb, numAwgnBpskBer);
hold on;
semilogy(snrDb, numZfBpskBer);
hold on;
semilogy(snrDb, numMrcBpskBer);
hold on;
semilogy(snrDb, numMrtBpskBer);
hold on;
semilogy(snrDb, numAlamoutiBpskBer);
grid on;
legend('AWGN', 'ZF', '2-MRC', '2-MRT', '2-Alamouti', 'location', 'sw');
xlabel('SNR per bit (dB)');
ylabel('BER');
