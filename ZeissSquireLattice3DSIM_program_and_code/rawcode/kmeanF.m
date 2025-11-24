function [kA,coe] =  kmeanF(obj)
npeak=(size(obj.Separated,3)-1)/2;
kA = zeros(npeak,2);
coe=0;
for i = 1:npeak
    fSa1Tnoisy0 = obj.Separated(:,:,1);
    fSa1Tnoisyp = obj.Separated(:,:,(i-1)*2+2);
    kA(i,:) = IlluminationFreqTIRF(obj,fSa1Tnoisy0,fSa1Tnoisyp,obj.kAo(i,:));
    if npeak>3
        if i ==1||i ==5%1级
            coe=coe+norm(kA(i,:))./npeak*2;
            % kaotest=norm(obj.kAo(i,:))./npeak*2
            katest=norm(kA(i,:))*2
        elseif i ==2||i ==3 %2级
            coe=coe+norm(kA(i,:))./npeak;
            % kaotest=norm(obj.kAo(i,:))./npeak
            katest=norm(kA(i,:))
        else%1.5级
            coe=coe+norm(kA(i,:))./npeak*sqrt(2);
            % kaotest=norm(obj.kAo(i,:))./npeak*sqrt(2)
            katest=norm(kA(i,:))*sqrt(2)
        end
    else
        coe=coe+norm(kA(i,:))./npeak;
    end
end
coe=coe/obj.Kotf;
end