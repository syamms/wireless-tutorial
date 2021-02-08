function [deSymbol] = comb_alamouti(rxSymbol, fading)
	% Function:
	%   - 1- or 2-rx space-time Alamouti receiver
	%
	% Input:
	%   - rxSymbol: received symbols over all receive antennas
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - deSymbol: filtered symbol for detection
	%
	% Comments:
	%	- constrained to 2-tx transmitter and 1- or 2-rx receiver
	%   - encode multiple symbol copies over different transmit antennas to combat fading
	%   - assume channel is unchanged over 2 consecutive slots
	%   - require CSIR, pair with Alamouti precoder at the transmitter
	%   - full receive array gain, no transmit array gain
	%   - fixed capacity, enhanced BER, 3dB worse than MRC/MRT
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	nRxs = size(fading, 1);


	switch nRxs
	case 1
		effFading = [fading; conj(fading(2)), - conj(fading(1))];
		effSymbol = reshape(rxSymbol, [2, length(rxSymbol) / 2]);
		effSymbol(2, :) = conj(effSymbol(2, :));
	case 2
		effFading = [fading; conj(fading(1, 2)), - conj(fading(1, 1)); conj(fading(2, 2)), - conj(fading(2, 1))];
		effSymbol = reshape(rxSymbol, [4, length(rxSymbol) / 2]);
		effSymbol(3 : 4, :) = conj(effSymbol(3 : 4, :));
	end
	deSymbol = reshape(effFading' * effSymbol, [1, length(rxSymbol)]);
end
