function [fading] = fading_kronecker(nTxs, nRxs, corTx, corRx)
	% Function:
	%   - Kronecker model for spatially correlated Rayleigh fading channels
	%
	% Input:
	%   - nTxs: number of tranmit antennas
	%	- nRxs: number of receive antennas
	%	- corTx: transmit correlation coefficient
	%	- corRx: receive correlation coefficient
	%
	% Output:
	%   - fading: spatially correlated Rayleigh fading
	%
	% Comments:
	%   - the correlation coefficients only influence channel strength distribution
	%
	% Author & Date: Yang (i@snowztail.com) - 16 Mar 19

	fading = toeplitz(corRx .^ (0 : nRxs - 1)) ^ (1 / 2) * sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs)) * toeplitz(corTx .^ (0 : nTxs - 1)) ^ (1 / 2);
end
