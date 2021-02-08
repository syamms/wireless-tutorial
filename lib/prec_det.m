function [txSymbol] = prec_det(symbol, fading)
	% Function:
	%   - dominant eigenmode transmission precoder
	%
	% Input:
	%   - symbol: symbol to be transmitted on multiple transmit antennas
	%   - fading: channel fading coefficient
	%
	% Output:
	%   - txSymbol: filtered symbol for transmission
	%
	% Comments:
	%   - decompose the MIMO channels into parallel links
	%   - allocate all power to the strongest link to transmit one data stream
	%   - require CSIT, pair with DET combiner at the receiver
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	[~, ~, v] = svd(fading);
	detPrecoder = v(:, 1);

	txSymbol = detPrecoder * symbol;
end
