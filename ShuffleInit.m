%%  ShuffleInit
%   Permutation of matrix generation
%      AB....YZ --> ZAB....Y
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: nothing
%   Copyright @ 2012-2017 Hermann Kampermann
%   Email: Hermann.Kampermann@hhu.de
%   Last updated: August 06, 2017

function iMap = ShuffleInit(nPart, dimHL)

iMap(1:dimHL^nPart) = 0;
stateSii(1:nPart) = ' ';
for ii = 1:dimHL^nPart
    stateii = dec2base(ii-1,dimHL,nPart);
    stateSii(1) = stateii(nPart);
    stateSii(2:nPart) = stateii(1:nPart-1);
    iMap(ii) = base2dec(stateSii,dimHL) + 1;
end

end

