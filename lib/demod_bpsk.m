function [bpskBit, bpskSymbolPower] = demod_bpsk(bpskSymbol)
	% Function:
	%   - maximum-likelihood detector for BPSK symbols
	%
	% Input:
	%   - bpskSymbol: BPSK symbol stream
	%
	% Output:
	%   - bpskBit: detected bit stream
	%   - bpskPower: average BPSK symbol power
	%
	% Comments:
	%   - signal space is actually 2-d due to complex noise
	%   - real part positive -> symbol closer to sqrt(E_s) -> bit 0
	%   - real part negative -> symbol closer to -sqrt(E_s) -> bit 1
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	bpskBit = 1 / 2 * (1 - sign(real(bpskSymbol)));
	bpskSymbolPower = norm(bpskSymbol) ^ 2 / length(bpskSymbol);
end
