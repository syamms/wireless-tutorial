function [deSymbol] = comb_det(rxSymbol, fading)
	% Function:
	%   - dominant eigenmode transmission combiner
	%
	% Input:
	%   - rxSymbol: received symbols over all receive antennas
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - deSymbol: filtered symbol for detection
	%
	% Comments:
	%   - decompose the MIMO channels into parallel links
	%   - allocate all power to the strongest link to transmit one data stream
	%   - require CSIR, pair with DET precoder at the transmitter
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	[u, ~, ~] = svd(fading);
	detCombiner = u(:, 1)';

	deSymbol = detCombiner * rxSymbol;
end
