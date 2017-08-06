%%  RandomLocalUnitary
%     Random local unitary (uniform w.r.t. the Haar measure)
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: nothing
%   Copyright @ 2012-2017 Hermann Kampermann
%   Email: Hermann.Kampermann@hhu.de
%   Last updated: August 06, 2017

function WPrimeRand = RandomLocalUnitary(WPrimeInit, nPart, dimHL)

[Q1,R] = qr(random('Uniform',0,1,dimHL,dimHL) + 1j*random('Uniform',0,1,dimHL,dimHL));
  Udum = Q1*diag(sign(diag(R)));
 
% qr --- Orthogonal-triangular decomposition
% [Q,R] = qr(A), where A is m-by-n, produces an m-by-n upper triangular
% matrix R and an m-by-m unitary matrix Q so that A = Q*R.
for ii = 1:nPart-1
   [Q1,R] = qr(random('Uniform',0,1,dimHL,dimHL) + 1j*random('Uniform',0,1,dimHL,dimHL));
       Q1 = Q1*diag(sign(diag(R)));
     Udum = kron(Q1,Udum);
end

WPrimeRand = Udum*WPrimeInit;

