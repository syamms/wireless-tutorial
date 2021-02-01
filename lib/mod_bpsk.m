function [bpskSymbol] = mod_bpsk(bitStream, symbolPower)
	% Function:
	%   - map bit stream to uncoded BPSK symbols
	%
	% Input:
	%   - bitStream: bit stream in 0 and 1
	%   - symbolPower: average symbol power
	%
	% Output:
	%   - bpskSymbol: uncoded BPSK symbols
	%
	% Restraints:
	%   - plain symbol without error correction coding
	%
	% Comments:
	%   - signal space is 1-D
	%   - assume initial phase is 0
	%   - 0 -> sqrt(E_s), 1 -> -sqrt(E_s)
	%
	% Author & Date: Yang (i@snowztail.com) - 21 Jan 19

	bpskSymbol = sqrt(symbolPower) * (1 - 2 * bitStream);
end
