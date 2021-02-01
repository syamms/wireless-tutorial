function [deSymbol] = comb_zf(rxSymbol, fading)
	% Function:
	%   - zero-forcing combining with single or multiple receive antennas
	%
	% Input:
	%   - rxSymbol: received symbols over all receive antennas
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - deSymbol: filtered symbol for detection
	%
	% Comments:
	%   - lead to noise amplification problem at low SNR
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	deSymbol = rxSymbol / fading;
end
