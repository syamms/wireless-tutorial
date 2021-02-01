function [qpskSymbol] = mod_qpsk(bitStream, symbolPower)
	% Function:
	%   - map bit stream to uncoded QPSK symbols
	%
	% Input:
	%   - bitStream: bit stream in 0 and 1
	%   - symbolPower: average symbol power
	%
	% Output:
	%   - qpskSymbol: uncoded QPSK symbols
	%
	% Restraints:
	%   - plain symbol without error correction coding
	%
	% Comments:
	%   - signal space is 2-D
	%   - assume initial phase is pi / 4 (symbols on the diagonal of quadrants)
	%   - [0, 0] -> [sqrt(E_s / 2), sqrt(E_s / 2)]
	%   - [0, 1] -> [sqrt(E_s / 2), -sqrt(E_s / 2)]
	%   - [1, 0] -> [-sqrt(E_s / 2), sqrt(E_s / 2)]
	%   - [1, 1] -> [-sqrt(E_s / 2), -sqrt(E_s / 2)]
	%
	% Author & Date: Yang (i@snowztail.com) - 21 Jan 19

	bitStream = 1 - 2 * reshape(bitStream, [2, length(bitStream) / 2]);
	qpskSymbol = sqrt(symbolPower / 2) * (bitStream(1, :) + 1i * bitStream(2, :));
end
