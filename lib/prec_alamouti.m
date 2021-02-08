function [txSymbol] = prec_alamouti(symbol)
	% Function:
	%   - 2-tx space-time Alamouti precoder
	%
	% Input:
	%   - symbol: symbol to be transmitted on multiple transmit antennas
	%
	% Output:
	%   - txSymbol: filtered symbol for transmission
	%
	% Comments:
	%	- constrained to 2-tx transmitter
	%   - encode multiple symbol copies over different transmit antennas to combat fading
	%   - assume channel is unchanged over 2 consecutive slots
	%   - does not require CSIT, pair with Alamouti combiner at the receiver
	%   - fixed capacity, enhanced BER, 3dB worse than MRC/MRT
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	txSymbol = zeros(2, length(symbol));

	txSymbol(1, 1 : 2 : end) = sqrt(1 / 2) * symbol(1 : 2 : end);
	txSymbol(1, 2 : 2 : end) = - sqrt(1 / 2) * symbol(2 : 2 : end)';

	txSymbol(2, 1 : 2 : end) = sqrt(1 / 2) * symbol(2 : 2 : end);
	txSymbol(2, 2 : 2 : end) = sqrt(1 / 2) * symbol(1 : 2 : end)';
end
