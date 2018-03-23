function[Z1,Z2] = improper(sz,rz,N)
% This function generates two complex-valued time series of length N as
% specified by sz, rz, and N where:
%   sz:     is the autocovariance from lag 0 to N (will be embedded with zeros
%           up to lag N if unspecified)
%   rz:     is the complementary covariance (will be embedded with zeros
%           up to lag N if unspecified)
%   N:      is the desired length of the sequence
a1=length(sz); a2=length(rz);
if min(a1,a2)<N+1
    disp('Warning: length of specified covariances are less than lag N, remaining lags up to N have been embedded with zeros')
    if a1<N+1
        sz(a1+1:N+1)=0;
    end
    if a2<N+1
        rz(a2+1:N+1)=0;
    end
end
% lines 2 to 5 in Algorithm 1:
sxx = real(0.5*(sz(1:N+1)+rz(1:N+1)));
syy = real(0.5*(sz(1:N+1)-rz(1:N+1)));
sxy = imag(0.5*(rz(1:N+1)-sz(1:N+1)));
syx = imag(0.5*(rz(1:N+1)+sz(1:N+1)));
% lines 6 to 8 in Algorithm 1:
cxx = [sxx sxx(N:-1:2)];
cyy = [syy syy(N:-1:2)];
cxy = [syx sxy(N:-1:2)];
% line 9 of Algorithm 1:
Lxx = real(fft(cxx)); Lyy = real(fft(cyy)); Lxy = fft(cxy);
if min(Lxx)< 0
    disp('Warning: negative eigenvalues in circulant matrix, these have been set to zero, method is now approximate (you can try with a larger N to see if this is exact)')
    Lxx = max(0,Lxx);
end
if min(Lyy)< 0
    disp('Warning: negative eigenvalues in circulant matrix, these have been set to zero, method is now approximate (you can try with a larger N to see if this is exact)')
    Lyy = max(0,Lyy);
end
% Line 10 of Algorithm 1:
WW=(1/sqrt(2))*(randn(2,2*N)+1i*randn(2,2*N));
% Lines 11 to 14 in Algorithm 1:
AA=ones(2,2); AA=2*AA;
ZZ = zeros(2,2*N);
for ii = 1:2*N;
    if Lxx(ii)>0 && Lyy(ii)>0
        AO=2*Lxy(ii)./sqrt(Lxx(ii).*Lyy(ii));
    else
        AO=0;
    end
    AA(1,2)=AO;
    AA(2,1)=conj(AO);
    [EVE,EVA]=eig(AA);
    if min(diag(EVA))<0
        disp('Warning: negative eigenvalues in Cholesky matrix, these have been set to zero, method is now approximate (you can try with a larger N to see if this is exact)')
        AA = EVE*max(EVA,0.000001)/EVE; % set eigenvalues to very close to zero such that matrix is positive definite
        AA(1,1)=abs(AA(1,1)); AA(2,2)=abs(AA(2,2));
    end
    SIG=chol(AA)';
    ZZ(:,ii)=SIG*WW(:,ii);
end
% Line 15 of Algorithm 1:
XX = fft(sqrt(Lxx).*ZZ(1,:))/sqrt(2*N);
YY = fft(sqrt(Lyy).*ZZ(2,:))/sqrt(2*N);
% Line 16 of Algorithm 1:
X1 = real(XX(1:N)); Y1 = real(YY(1:N));
Z1=X1+1i.*Y1;
% Line 17 of Algorithm 1:
X2 = imag(XX(1:N)); Y2 = imag(YY(1:N));
Z2=X2+1i.*Y2;