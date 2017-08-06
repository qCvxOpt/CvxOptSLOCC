%%  WSLOCCMaxG
%      W-state SLOCC maximization (n-parties)
%
%   URL: https://github.com/qCvxOpt/CvxOptSLOCC
%   Requires: nothing
%   Copyright @ 2012-2017 Hermann Kampermann
%   Email: Hermann.Kampermann@hhu.de
%   Last updated: August 06, 2017

function WPrime = WSLOCCMaxG(Rho, WPrimeInit, iMap, dimHL, nPart) 

iPart = 50;
rhos = Rho;
WPrime = WPrimeInit;
dHL2 = dimHL*dimHL;
dges = dimHL^nPart;
dME = dimHL^(nPart-1);
ed = 1.;
lp = 0;
while ((ed>10^-14||iP<nPart) && lp<iPart*nPart)

iP = mod(lp,nPart) + 1;
lp = lp + 1;

% Here we need to do the state reshuffeling routine (exchange parties, 
% such that the considered party is at the first place
% if iP==nPart we have the original ordering
% AB...YZ --> ZAB...Y
[WPrime,rhos] = ShuffleR(WPrime,rhos,iMap,dges);

Cnorm(1:dHL2,1:dHL2) = zeros;
for aa = 1:dimHL
    for bb = 1:dimHL
        for cc = 1:dME   
            aac = (aa-1)*dME+cc;
            bbc = (bb-1)*dME+cc;
            Cnorm(aa,bb) = Cnorm(aa,bb) + conj(WPrime(aac))*WPrime(bbc);                                                                            
        end
    end
end

% Here we copy the blocks
for ii = 1:dimHL-1
    for aa = 1:dimHL
        for bb = 1:dimHL
            Cnorm(ii*dimHL+aa,ii*dimHL+bb) = Cnorm(aa,bb);
        end 
    end
end

%
WRW(1:dHL2,1:dHL2) = zeros;
for ii = 1:dimHL
    for jj = 1:dimHL
        for kk = 1:dimHL
            for ll = 1:dimHL
                aa = (ii-1)*dimHL+jj;
                bb = (kk-1)*dimHL+ll;
                for cc = 1:dME   
                    for dd = 1:dME
                        wjc = (jj-1)*dME+cc;
                        wld = (ll-1)*dME+dd;
                        ric = (ii-1)*dME+cc;
                        rkd = (kk-1)*dME+dd;
                        WRW(aa,bb) = WRW(aa,bb)+...
                                conj(WPrime(wjc))*rhos(ric,rkd)*WPrime(wld);
                    end
                end
            end    
        end
    end
end

[U,E] = eig(WRW-real(WPrime'*rhos*WPrime)*Cnorm);
[ed,imm] = max(real(diag(E)));

if (ed>10^-14)
  a(1:dimHL,1:dimHL) = zeros;
  iic = 0;
  for ii = 1:dimHL
      for jj = 1:dimHL
          iic = iic+1;
          a(ii,jj) = U(iic,imm);
      end 
  end
  a = a + 3.*ed*eye(dimHL);
  
  WPrime = kron(a,eye(dME))*WPrime;
  nor = sqrt(real(WPrime'*WPrime));
  WPrime = WPrime/nor; 
  
% !!!
  if (nor<10^-11)
      'Problem in SLOCC'
      nor,lp,iP,ed
      a
      mat = real(diag(E));
      mat
      pause
  end
% !!!
end

end

end


%% Permutation of the qudit indices
%     AB....YZ --> ZAB....Y
function [WPrime, rhos] = ShuffleR(WPrime, rhos, iMap, dges)

Wd = WPrime;
Rd = rhos;
for ii = 1:dges
    WPrime(iMap(ii)) = Wd(ii);
    for jj = 1:dges
        rhos(iMap(ii),iMap(jj)) = Rd(ii,jj);
    end
end

end

