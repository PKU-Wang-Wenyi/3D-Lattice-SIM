function [separate]= PCMseparate_tempF_3D(obj)
fS1aTnoisy = obj.fSnoisy(:,:,1);
fS2aTnoisy = obj.fSnoisy(:,:,2);
fS3aTnoisy = obj.fSnoisy(:,:,3);
fS4aTnoisy = obj.fSnoisy(:,:,4);
fS5aTnoisy = obj.fSnoisy(:,:,5);
fS6aTnoisy = obj.fSnoisy(:,:,6);
fS7aTnoisy = obj.fSnoisy(:,:,7);
fS8aTnoisy = obj.fSnoisy(:,:,8);
fS9aTnoisy = obj.fSnoisy(:,:,9);
fS10aTnoisy = obj.fSnoisy(:,:,10);
fS11aTnoisy = obj.fSnoisy(:,:,11);
fS12aTnoisy = obj.fSnoisy(:,:,12);
fS13aTnoisy = obj.fSnoisy(:,:,13);
fDuplex2 = fS2aTnoisy - fS1aTnoisy;
fDuplex3 = fS3aTnoisy - fS1aTnoisy;
fDuplex4 = fS4aTnoisy - fS1aTnoisy;
fDuplex5 = fS5aTnoisy - fS1aTnoisy;
fDuplex6 = fS6aTnoisy - fS1aTnoisy;
fDuplex7 = fS7aTnoisy - fS1aTnoisy;
fDuplex8 = fS8aTnoisy - fS1aTnoisy;
fDuplex9 = fS9aTnoisy - fS1aTnoisy;
fDuplex10 = fS10aTnoisy - fS1aTnoisy;
fDuplex11 = fS11aTnoisy - fS1aTnoisy;
fDuplex12 = fS12aTnoisy - fS1aTnoisy;
fDuplex13 = fS13aTnoisy - fS1aTnoisy;
Kai2Opt3D0 = @(phase0)Kai2Opt3D(obj,phase0,fDuplex2,fDuplex3,fDuplex4,fDuplex5,fDuplex6,fDuplex7,fDuplex8,fDuplex9,fDuplex10,fDuplex11,fDuplex12,fDuplex13);
options = optimset('LargeScale','off','Algorithm',...
    'active-set','MaxFunEvals',4000,'MaxIter',2000,'Display','notify');
phase0 = gpuArray([2*pi/13 4*pi/13 6*pi/13 8*pi/13 10*pi/13 12*pi/13 14*pi/13 16*pi/13 18*pi/13 20*pi/13 22*pi/13 24*pi/13]); % initial guess
[phaseShift,fval] = fminsearch(Kai2Opt3D0,phase0,options);
phaseShift_ini = 0;
[fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2,fDp3,fDm3,fDp4,fDm4,fDp5,fDm5] = SeparatedComponents3D(...
    phaseShift,phaseShift_ini,fS1aTnoisy,fS2aTnoisy,fS3aTnoisy,fS4aTnoisy,fS5aTnoisy,fS6aTnoisy,fS7aTnoisy,fS8aTnoisy,fS9aTnoisy,fS10aTnoisy,fS11aTnoisy,fS12aTnoisy,fS13aTnoisy);
separate=cat(3,fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2,fDp3,fDm3,fDp4,fDm4,fDp5,fDm5);
% figure;mesh(log(abs(fDo)+1))
% figure;mesh(log(abs(fDp)+1))
end