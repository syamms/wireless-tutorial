function [deSymbol] = comb_mrc(rxSymbol, fading)
	% Function:
	%   - maximal ratio combining with multiple receive antennas
	%
	% Input:
	%   - rxSymbol: received symbols over all receive antennas
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - deSymbol: filtered symbol for detection
	%
	% Comments:
	%   - enhance array and diversity gains, require CSIR
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	deSymbol = fading' / norm(fading) * rxSymbol;
end
