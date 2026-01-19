function [obj,coe] = KapproxEstimationF_PCA(obj,param)
numpixelsx=obj.w;
cnt=[numpixelsx/2+1,numpixelsx/2+1];
cutoff=1000/(0.5*param.lambda/param.NA);
cyclesPerMicron = (1./(numpixelsx*param.pixelsize*0.001));
[x,y]=meshgrid(1:numpixelsx,1:numpixelsx);
rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);
Mask1=double(rad<=1.0*(cutoff/cyclesPerMicron+1));
% NotchFilter0=getotfAtt(numpixelsx,cyclesPerMicron,0.5*cutoff,0,0);%这里必须精细的调节mask，要不然就会跑偏
NotchFilter0=getotfAtt(numpixelsx,cyclesPerMicron,0.3*cutoff,0,0);
NotchFilter1=NotchFilter0.*Mask1;
Mask2=double(rad<=1.10*(cutoff/cyclesPerMicron+1));
NotchFilter2=NotchFilter0.*Mask2;
%%
numangles=6;
nrBands = 13;
phaOff=0;
fac = ones(1,nrBands);
separateII = separateBandshifi_lattice(obj.fSnoisy,phaOff,nrBands,fac);
%%
CrossCorrelation=zeros(size(Mask2,1),size(Mask2,2),numangles);
k0=zeros(1,numangles);
for I=1:numangles
    c0 = separateII(:,:,1);
    c3 = separateII(:,:,2*I);
    c0 = c0./(max(max(abs(separateII(:,:,1))))).*NotchFilter1;
    c3 = c3./(max(max(abs(separateII(:,:,2*I))))).*NotchFilter1;
    c0 = FFT2D(c0,false);
    c3 = FFT2D(c3,false);
    c3 = c3.*conj(c0);
    c3 = c3./max(max(c3));
    vec = fftshift(FFT2D(c3,true));
    CrossCorrelation(:,:,I) = vec;
    temp = vec.*NotchFilter2;
    temp = log(1+abs(temp));
    temp = temp./max(max(temp));
    % figure;mesh(temp)
    [yPos,xPos] = find(temp==max(max(temp)));
    peak.xPos = xPos(1);
    peak.yPos = yPos(1);
    k0(I) = sqrt((peak.xPos-cnt(1))^2+(peak.yPos-cnt(2))^2);
end
coe=mean(k0)/obj.Kotf;
%% PCA
Filter_size=11;
Mask_size=3;
PCA_num=1;
overlap=0.15;
for I = 1:numangles
    vec = CrossCorrelation(:,:,I);
    temp = vec.*NotchFilter2;
    temp = log(1+abs(temp));
    temp = temp./max(max(temp));

    [yPos,xPos] = find(temp==max(max(temp)));
    peak.xPos = xPos(1);
    peak.yPos = yPos(1);
    for pca_size=1:PCA_num
        if(pca_size==1)
            kx = (peak.xPos-cnt(2)) ;
            ky = (peak.yPos-cnt(1));
            old_kx=kx ;
            old_ky=ky;
        else
            kx = old_kx ;
            ky = old_ky;
        end
        Temp=NfourierShift(placeFreq(separateII(:,:,2*I)),-(2-1)*old_kx,-(2-1)* old_ky);
        MASK = ones(size(Temp));
        MASK(numpixelsx-Mask_size+1:numpixelsx+Mask_size+1, numpixelsx-Mask_size+1:numpixelsx +Mask_size+1) = 0;
        Temp(MASK==1) = 0;
        ROI = Temp(numpixelsx-Filter_size+1:numpixelsx+Filter_size+1, numpixelsx-Filter_size+1:numpixelsx+Filter_size+1);
        SIZE = 2*Filter_size+1;
        Space_ROI=(ifft2(ifftshift(ROI)));
        Phase_ROI=exp(1i*angle(Space_ROI));                       %加掩码算子划定ROI

        [U,S,V] = svd(Phase_ROI);                                                  % 奇异值分解
        Y = U*S*V';                                                                               %主成分分布
        SS=zeros(SIZE,SIZE);
        SS(1,1)=S(1,1);                                                                         %提取第一主成分
        Phase_ROI= U*SS*V';

        UNV1=unwrap(angle(V(:,1)));                                                % 对右奇异矩阵V的第一列元素进行相位展开，通过最小二乘法进行线性拟合，拟合出的斜率即为残余的亚像素波矢ksub在竖直方向上的分量；
        UNV_T=UNV1';
        UNV_T_fit=UNV_T(1,2:SIZE);

        FIT=polyfit(1:SIZE-1,UNV_T_fit,1);
        dx = - FIT(1)*SIZE/2/pi;                                                              %波矢ksub在竖直方向上的分量dx
        NewUNV1=polyval(FIT,1:SIZE);
        NEWV1=exp(1i.*((NewUNV1')));
        V(:,1)=NEWV1;

        UNU1=unwrap(angle(U(:,1)));                                                %对左奇异矩阵U的第一行元素进行相位展开，通过最小二乘法进行线性拟合，拟合出的斜率即为残余的亚像素波矢ksub在水平方向上的分量
        UNN_T=UNU1';
        UNN_T_fit=UNN_T(1,2:SIZE);
        FIT=polyfit(1:SIZE-1,UNN_T_fit,1);
        dy = FIT(1)*SIZE/2/pi;                                                                 %波矢ksub在水平方向上的分量dy
        NewUN1=polyval(FIT,1:SIZE);
        NEWU1=exp(1i.*(NewUN1'));
        U(:,1)=NEWU1;

        NEW_Phase_ROI= U*SS*V';

        old_kx= -(old_kx + dx) ;
        old_ky= -(old_ky + dy);

        % p1 = getPeak(separateII(:,:,1),separateII(:,:,lb),0,1,OtfProvider,old_kx/2,old_ky/2,overlap);
        p = getPeak(separateII(:,:,1)/(max(max(abs(separateII(:,:,1))))),separateII(:,:,2*I)/(max(max(abs(separateII(:,:,1))))),0,2,param.OtfProvider1,old_kx,old_ky,overlap);

        params.Dir(I).px = -old_kx/2;
        params.Dir(I).py = -old_ky/2;
        obj.kAmean(I,:)=[old_ky,old_kx];
        obj.phaseA(I)= - angle(NEW_Phase_ROI(Filter_size+1,Filter_size+1));
        Temp_m1 = abs(p);
        obj.modFac(I) = Temp_m1;

    end
end
obj.phaseA=obj.phaseA';
obj.modFac=obj.modFac';
end
function [out] = placeFreq( in )                                          %频谱搬移
siz=size(in);
w=siz(2);
h=siz(1);
out=zeros(2*w,2*h,'like',in);
out(h/2+1:h+h/2,w/2+1:w+w/2)=in;
end

function [outv] = NfourierShift( inv,kx,ky )
inv=(ifft2(fftshift(inv)));
siz=size(inv);
[x,y]=meshgrid(0:siz(2)-1,0:siz(1)-1);
x=x/siz(2);
y=y/siz(1);
outv=inv.*exp(2*pi*1i*(ky*y+kx*x));
outv=fftshift(fft2((outv)));
end
