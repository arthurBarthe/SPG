function [out]= lkh_(x, SX, ker, delta,LB,UB,MF,ZEROF)
%This function computes the likelihood defined in the paper, using the
%kernel sent as the argument called ker. The latter needs to be precomputed
%via the call of the function kernel. 
%
%Note: the coriolis frequencies are accounted for in the kernel k (which is
%complex-valued). But we still use 3 args when calling the function
%complexouacvs. The second arg is no longer the peak of the OU but a shift
%from the coriolis frequencies. For that reason the initial guess xb2(2) is
%to be set close to zero (no shift from the coriolis frequencies).
% acv = maternacvs(x(4:6).*xb(4:6),N,delta) .* (1-(0:N-1)/N);
% acv = acv + complexouacvs(x(1:3) .* xb(1:3), N, delta) .* ker;
% ESF2=2*fft(acv)-acv(1); ESF3=abs(real(fftshift(ESF2)));
ESF3 = S_(ker, x, delta);
if ZEROF==0
    out=sum(log(ESF3([LB:MF-1 MF+1:UB])))+sum(SX([LB:MF-1 MF+1:UB])./ESF3([LB:MF-1 MF+1:UB]));
else
    out=sum(log(ESF3([LB:UB])))+sum(SX([LB:UB])./ESF3([LB:UB]));
end
end
