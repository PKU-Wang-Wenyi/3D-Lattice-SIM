function obj=PCMseparateF_3D(obj)
phase = obj.phaseA;
%% Separating the three frequency components
MF = 1.0;
%% Transformation Matrix
M = 0.5*[1 0.5*MF*exp(-1i*phase(1)) 0.5*MF*exp(+1i*phase(1)) 0.5*MF*exp(-1i*phase(2)) 0.5*MF*exp(+1i*phase(2)) 0.5*MF*exp(-1i*phase(3)) 0.5*MF*exp(+1i*phase(3)) 0.5*MF*exp(-1i*phase(4)) 0.5*MF*exp(+1i*phase(4)) 0.5*MF*exp(-1i*phase(5)) 0.5*MF*exp(+1i*phase(5)) 0.5*MF*exp(-1i*phase(6)) 0.5*MF*exp(+1i*phase(6));
    1 0.5*MF*exp(-1i*(phase(1)+2*pi/13)) 0.5*MF*exp(+1i*(phase(1)+2*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*2*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*2*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*2*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*2*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*2*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*2*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*2*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*2*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*2*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*2*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+4*pi/13)) 0.5*MF*exp(+1i*(phase(1)+4*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*4*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*4*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*4*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*4*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*4*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*4*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*4*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*4*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*4*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*4*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+6*pi/13)) 0.5*MF*exp(+1i*(phase(1)+6*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*6*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*6*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*6*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*6*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*6*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*6*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*6*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*6*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*6*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*6*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+8*pi/13)) 0.5*MF*exp(+1i*(phase(1)+8*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*8*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*8*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*8*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*8*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*8*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*8*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*8*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*8*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*8*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*8*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+10*pi/13)) 0.5*MF*exp(+1i*(phase(1)+10*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*10*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*10*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*10*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*10*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*10*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*10*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*10*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*10*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*10*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*10*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+12*pi/13)) 0.5*MF*exp(+1i*(phase(1)+12*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*12*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*12*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*12*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*12*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*12*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*12*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*12*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*12*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*12*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*12*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+14*pi/13)) 0.5*MF*exp(+1i*(phase(1)+14*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*14*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*14*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*14*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*14*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*14*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*14*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*14*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*14*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*14*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*14*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+16*pi/13)) 0.5*MF*exp(+1i*(phase(1)+16*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*16*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*16*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*16*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*16*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*16*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*16*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*16*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*16*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*16*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*16*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+18*pi/13)) 0.5*MF*exp(+1i*(phase(1)+18*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*18*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*18*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*18*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*18*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*18*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*18*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*18*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*18*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*18*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*18*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+20*pi/13)) 0.5*MF*exp(+1i*(phase(1)+20*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*20*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*20*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*20*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*20*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*20*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*20*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*20*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*20*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*20*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*20*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+22*pi/13)) 0.5*MF*exp(+1i*(phase(1)+22*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*22*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*22*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*22*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*22*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*22*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*22*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*22*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*22*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*22*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*22*pi/13));
    1 0.5*MF*exp(-1i*(phase(1)+24*pi/13)) 0.5*MF*exp(+1i*(phase(1)+24*pi/13)) 0.5*MF*exp(-1i*(phase(2)+2*24*pi/13)) 0.5*MF*exp(+1i*(phase(2)+2*24*pi/13)) 0.5*MF*exp(-1i*(phase(3)+3*24*pi/13)) 0.5*MF*exp(+1i*(phase(3)+3*24*pi/13)) 0.5*MF*exp(-1i*(phase(4)+4*24*pi/13)) 0.5*MF*exp(+1i*(phase(4)+4*24*pi/13)) 0.5*MF*exp(-1i*(phase(5)+5*24*pi/13)) 0.5*MF*exp(+1i*(phase(5)+5*24*pi/13)) 0.5*MF*exp(-1i*(phase(6)+6*24*pi/13)) 0.5*MF*exp(+1i*(phase(6)+6*24*pi/13));];
Minv = pinv(M);
%%
[H, W, ~] = size(obj.fSnoisy(:,:,1));
fSeparated = zeros(H, W, 13, 'like', obj.fSnoisy);

for k = 1:13
    for j = 1:13
        fSeparated(:, :, k) = fSeparated(:, :, k) + Minv(k, j) * obj.fSnoisy(:, :, j);
    end
end

obj.Separated = fSeparated;

% figure;mesh(log(abs(fSeparated(:,:,1)+1)))
% figure;mesh(log(abs(fSeparated(:,:,2)+1)))
% figure;mesh(log(abs(fSeparated(:,:,4)+1)))
% figure;mesh(log(abs(fSeparated(:,:,6)+1)))
% figure;mesh(log(abs(fSeparated(:,:,8)+1)))
% figure;mesh(log(abs(fSeparated(:,:,10)+1)))
% figure;mesh(log(abs(fSeparated(:,:,12)+1)))



[H, W, ~, N] = size(obj.fSnoisy_all);
obj.Separated_all = zeros(H, W, 13, N, 'like', obj.fSnoisy_all);

for i = 1:N
    fSeparated = zeros(H, W, 13, 'like', obj.fSnoisy_all);
    for k = 1:13
        for j = 1:13
            fSeparated(:, :, k) = fSeparated(:, :, k) + Minv(k, j) * obj.fSnoisy_all(:, :, j, i);
        end
    end
    obj.Separated_all(:, :, :, i) = fSeparated;
end
% obj.fSnoisy=[];
% obj.fSnoisy_all=[];
% obj = rmfield(obj, {'fSnoisy', 'fSnoisy_all'});
end