function [qpskBit, qpskSymbolPower] = demod_qpsk(qpskSymbol)
	% Function:
	%   - maximum-likelihood detector for QPSK symbols
	%
	% Input:
	%   - qpskSymbol: QPSK symbol stream
	%
	% Output:
	%   - qpskBit: detected bit stream
	%   - qpskSymbolPower: average QPSK symbol power
	%
	% Comments:
	%   - signal space is 2-d
	%   - positive real part -> 1st bit closer to sqrt(E_s / 2) -> 1st bit 0
	%   - negative real part -> 1st bit closer to -sqrt(E_s / 2) -> 1st bit 1
	%   - positive imaginary part -> 2nd bit closer to i * sqrt(E_s / 2) -> 2nd bit 0
	%   - negative imaginary part -> 2nd bit closer to i * -sqrt(E_s / 2) -> 2nd bit 1
	%
	% Author & Date: Yang (i@snowztail.com) - 22 Jan 19

	qpskBit(1 : 2 : 2 * length(qpskSymbol) - 1) = 1 / 2 * (1 - sign(real(qpskSymbol)));
	qpskBit(2 : 2 : 2 * length(qpskSymbol)) = 1 / 2 * (1 - sign(imag(qpskSymbol)));
	qpskSymbolPower = norm(qpskSymbol) ^ 2 / length(qpskSymbol);
end
