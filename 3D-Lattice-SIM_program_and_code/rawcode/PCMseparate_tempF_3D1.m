function [separate] = PCMseparate_tempF_3D1(obj)
% === 提取所有傅里叶图像 ===
fS = obj.fSnoisy;  % [H, W, 13]
fS1 = fS(:,:,1);

% === 计算相对差分 fDuplex ===
fDuplex = fS(:,:,2:13) - fS1;  % [H, W, 12]
% === 定义优化目标函数 ===
G = obj.G; % [H, W]
Kai2Opt3D0 = @(phase0) Kai2Opt3D1(fDuplex, phase0,G);

% === 优化参数 ===
options = optimset('LargeScale','off', 'Algorithm','active-set', ...
                   'MaxFunEvals',4000, 'MaxIter',2000, 'Display','notify');

phase0 = gpuArray(2*pi*(1:12)/13); % 初始相位猜测
[phaseShift, fval] = fminsearch(Kai2Opt3D0, phase0, options);

% === 调用分离函数 ===
phaseShift_ini = 0;
% [fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2,fDp3,fDm3,fDp4,fDm4,fDp5,fDm5] = ...
%     SeparatedComponents3D(phaseShift, phaseShift_ini, fS(:,:,1:13));
[fS1aTnoisy, fS2aTnoisy, fS3aTnoisy, fS4aTnoisy, fS5aTnoisy, ...
 fS6aTnoisy, fS7aTnoisy, fS8aTnoisy, fS9aTnoisy, fS10aTnoisy, ...
 fS11aTnoisy, fS12aTnoisy, fS13aTnoisy] = deal(...
    fS(:,:,1), fS(:,:,2), fS(:,:,3), fS(:,:,4), fS(:,:,5), ...
    fS(:,:,6), fS(:,:,7), fS(:,:,8), fS(:,:,9), fS(:,:,10), ...
    fS(:,:,11), fS(:,:,12), fS(:,:,13));
[fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2,fDp3,fDm3,fDp4,fDm4,fDp5,fDm5] = SeparatedComponents3D(...
    phaseShift,phaseShift_ini,fS1aTnoisy,fS2aTnoisy,fS3aTnoisy,fS4aTnoisy,fS5aTnoisy,fS6aTnoisy,fS7aTnoisy,fS8aTnoisy,fS9aTnoisy,fS10aTnoisy,fS11aTnoisy,fS12aTnoisy,fS13aTnoisy);
separate = cat(3, fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2,fDp3,fDm3,fDp4,fDm4,fDp5,fDm5);

end
