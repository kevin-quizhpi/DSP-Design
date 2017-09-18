% qwrap.m - vectorized version of qwrap.c
%
% Usage: q_out = qwrap(D, q) = mod-(D+1) version of q
%
% D = maximum order, buffer length = D+1
% q = circular index into buffer, must be in the range -1<=q<=2*D
%
% q_out = wrapped version of q
%
% easier by anonymous function: Q = @(D,q) q + (D+1)*((q<0) - (q>D));

% S. J. Orfanidis, 332:447 DSPD, Fall 2017

function q_out = qwrap(D, q)

if nargin==0, help qwrap; return; end

q_out = q + (D+1)*((q<0) - (q>D));
