function [separate]= PCMseparate_tempF_2D(obj)
fS1aTnoisy = obj.fSnoisy(:,:,1);
fS2aTnoisy = obj.fSnoisy(:,:,2);
fS3aTnoisy = obj.fSnoisy(:,:,3);
fS4aTnoisy = obj.fSnoisy(:,:,4);
fS5aTnoisy = obj.fSnoisy(:,:,5);
fS6aTnoisy = obj.fSnoisy(:,:,6);
fS7aTnoisy = obj.fSnoisy(:,:,7);

fDuplex2 = fS2aTnoisy - fS1aTnoisy;
fDuplex3 = fS3aTnoisy - fS1aTnoisy;
fDuplex4 = fS4aTnoisy - fS1aTnoisy;
fDuplex5 = fS5aTnoisy - fS1aTnoisy;
fDuplex6 = fS6aTnoisy - fS1aTnoisy;
fDuplex7 = fS7aTnoisy - fS1aTnoisy;
Kai2Opt0 = @(phase0)Kai2Opt(obj,phase0,fDuplex2,fDuplex3,fDuplex4,fDuplex5,fDuplex6,fDuplex7);
options = optimset('LargeScale','off','Algorithm',...
    'active-set','MaxFunEvals',2000,'MaxIter',2000,'Display','notify');
phase0 = [2*pi/7 4*pi/7 6*pi/7 8*pi/7 10*pi/7 12*pi/7]; % initial guess
[phaseShift,fval] = fminsearch(Kai2Opt0,phase0,options);
phaseShift_ini = 0;
[fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2] = SeparatedComponents2D(...
    phaseShift,phaseShift_ini,fS1aTnoisy,fS2aTnoisy,fS3aTnoisy,fS4aTnoisy,fS5aTnoisy,fS6aTnoisy,fS7aTnoisy);
separate=cat(3,fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2);
% figure;mesh(log(abs(separate(:,:,1))+1))
% figure;mesh(log(abs(separate(:,:,2))+1))
% figure;mesh(log(abs(separate(:,:,4))+1))
% figure;mesh(log(abs(separate(:,:,6))+1))
end