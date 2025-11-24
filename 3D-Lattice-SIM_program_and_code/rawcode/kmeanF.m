function [kA,coe] =  kmeanF(obj)
npeak=(size(obj.Separated,3)-1)/2;
kA = zeros(npeak,2);
coe=0;
for i = 1:npeak
    fSa1Tnoisy0 = obj.Separated(:,:,1);
    fSa1Tnoisyp = obj.Separated(:,:,(i-1)*2+2);
    kA(i,:) = IlluminationFreqTIRF(obj,fSa1Tnoisy0,fSa1Tnoisyp,obj.kAo(i,:));
    if npeak>3
        if i ==1||i ==3||i ==4%1级
            coe=coe+norm(kA(i,:))./npeak*sqrt(3);
            % kaotest=norm(obj.kAo(i,:))./npeak*sqrt(3)
            katest=norm(kA(i,:))./npeak*sqrt(3)
        else%2级
            coe=coe+norm(kA(i,:))./npeak;
            % kaotest=norm(obj.kAo(i,:))./npeak
            katest=norm(kA(i,:))./npeak
        end
    else
        coe=coe+norm(kA(i,:))./npeak;
    end
end
coe=coe/obj.Kotf;
end