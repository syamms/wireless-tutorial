%% * BPSK/QPSK over AWGN channel
clear; close all; clc; setup;
nTxs = 1;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 10;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e6;
bitStream = round(rand(1, nBits));

[bpskSymbol] = mod_bpsk(bitStream, bitPower);
[qpskSymbol] = mod_qpsk(bitStream, 2 * bitPower);

anaBer = qfunc(sqrt(2 * db2pow(snrDb)));
numBpskBer = zeros(length(snrDb), 1);
numQpskBer = zeros(length(snrDb), 1);

for iSnr = 1 : length(snrDb)
	rxBpskSymbol = bpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
	rxQpskSymbol = qpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
	[bpskBit] = demod_bpsk(rxBpskSymbol);
	[qpskBit] = demod_qpsk(rxQpskSymbol);
	numBpskBer(iSnr) = sum(xor(bitStream, bpskBit)) / nBits;
	numQpskBer(iSnr) = sum(xor(bitStream, qpskBit)) / nBits;
end

figure('name', 'BER of BPSK/QPSK over an AWGN channel');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK');
xlabel('SNR per bit (dB)');
ylabel('BER');

%% * BPSK/QPSK over SISO Rayleigh fading channel
clear; close all; clc; setup;
nTxs = 1;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e4;
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
		rxBpskSymbol = fading * bpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		rxQpskSymbol = fading * qpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
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
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK');
xlabel('SNR per bit (dB)');
ylabel('BER');

%% * BPSK/QPSK over SIMO i.i.d. Rayleigh fading channel with 2-rx MRC combining
clear; close all; clc; setup;
nTxs = 1;
nRxs = 2;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e4;
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
		rxBpskSymbol = fading * bpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(bpskSymbol)) + 1i * randn(nRxs, length(bpskSymbol)));
		rxQpskSymbol = fading * qpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(qpskSymbol)) + 1i * randn(nRxs, length(qpskSymbol)));
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

figure('name', 'BER of BPSK/QPSK over SISO Rayleigh fading channel');
semilogy(snrDb, anaBer);
hold on;
semilogy(snrDb, anaBerUb);
hold on;
semilogy(snrDb, numBpskBer);
hold on;
semilogy(snrDb, numQpskBer);
grid on;
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK');
xlabel('SNR per bit (dB)');
ylabel('BER');

%% * BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx MRT beamforming
clear; close all; clc; setup;
nTxs = 2;
nRxs = 1;
bitPower = 1;
snrDb = 0 : 20;
noisePower = bitPower ./ db2pow(snrDb);
nBits = 1e3;
nChannels = 1e4;
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
		[txBpskSymbol] = tran_mrt(bpskSymbol, fading);
		[txQpskSymbol] = tran_mrt(qpskSymbol, fading);
		rxBpskSymbol = fading * txBpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(txBpskSymbol)) + 1i * randn(nRxs, length(txBpskSymbol)));
		rxQpskSymbol = fading * txQpskSymbol + sqrt(noisePower(iSnr) / 2) * (randn(nRxs, length(txQpskSymbol)) + 1i * randn(nRxs, length(txQpskSymbol)));
		[bpskBit] = demod_bpsk(rxBpskSymbol);
		[qpskBit] = demod_qpsk(rxQpskSymbol);
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
legend('Analytical: BPSK / QPSK', 'Upper Bound: BPSK / QPSK', 'Numerical: BPSK', 'Numerical: QPSK');
xlabel('SNR per bit (dB)');
ylabel('BER');

%% * Alamouti scheme with BPSK/QPSK over MISO i.i.d. Rayleigh fading channel with 2-tx
