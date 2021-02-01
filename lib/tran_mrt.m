function [txSymbol] = tran_mrt(symbol, fading)
	% Function:
	%   - maximal ratio transmission with multiple transmit antennas
	%
	% Input:
	%   - symbol: modulated symbols to be transmitted
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - txSymbol: filtered symbol for transmission
	%
	% Comments:
	%   - enhance array and diversity gains, require CSIT (could be challenging)
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	txSymbol = fading' / norm(fading) * symbol;
end
