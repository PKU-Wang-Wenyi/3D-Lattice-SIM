function [OTF_orders,OTF_every_slide] = get_modelOTF_fast(params)

% raw funtcion copyright Sjoerd Stallinga, TU Delft, 2017-2020.
% Wenyi Wang rewrite some function into GPU to save some time.  2023
%% Basic parameters
Nx = params.w;
Ny = params.w;
Nz = params.nz;

%% support size and sampling for OTF in frequency domain,
OTFparams.supportsizex = 1/2/params.pixelsize;
OTFparams.supportsizey = 1/2/params.pixelsize;
OTFparams.supportsizez = 1/2/params.Zstep;
OTFparams.Nsupportx = Nx;
OTFparams.Nsupporty = Ny;
OTFparams.Nsupportz = Nz;
OTFparams.shiftsupport = [1,1,0];

%% define struct OTFparams needed for the vectorial PSF model functions
OTFparams.refmed = params.medium_refractive_index;    % refractive index medium/specimen
OTFparams.refcov = params.coverglass_refractive_index;    % refractive index conver slip
OTFparams.refimm = params.Sample_refractive_index;    % refractive index immersion medium
OTFparams.refimmnom = OTFparams.refcov; % nominal value of refractive indexe immersion medium, for which the objective lens is designed
OTFparams.fwd = params.obj_workingdistence;          % free working distance from objective to cover slip
OTFparams.depth = 0;      % depth of imaged slice from the cover slip
OTFparams.NA = params.NA;            % NA of the objective lens
OTFparams.xemit = 0.0;                   % x-position focal point
OTFparams.yemit = 0.0;                   % y-position focal point
OTFparams.zemit = 0.0;                   % z-position focal point
OTFparams.ztype = 'medium';              % z-distances measured inside the medium

%% sampling real and spatial frequency space
OTFparams.Npupil = round(sqrt(Nx*Ny)); % sampling points in pupil plane
OTFparams.Mx = Nx;       % sampling points in image space in x
OTFparams.My = Ny;       % sampling points in image space in y
OTFparams.Mz = Nz;       % sampling points in image space in z (must be 1 for 2D)
OTFparams.xrange = Nx*params.pixelsize/2;           % 1/2-size of image space in x
OTFparams.yrange = Ny*params.pixelsize/2;           % 1/2-size of image space in x
OTFparams.zrange = [-Nz*params.Zstep/2,Nz*params.Zstep/2]; % [zmin,zmax] defines image space in z
OTFparams.pixelsize = params.pixelsize;        % pixel size in image space
OTFparams.samplingdistance = params.pixelsize; % sampling distance in image space
OTFparams.aberrations = [1,1,0.0; 1,-1,-0.0; 2,0,-0.0; 4,0,0.0; 2,-2,0.0; 2,2,0.0; 4,-2,0.0];
OTFparams.dipoletype = 'free'; % averaging over dipole orientations

%% loop over color channels
OTFparams.lambdaex = params.exlambda;                    % excitation wavelength
OTFparams.lambda = params.lambda;                      % emmission wavelength
OTFparams.aberrations(:,3) =  OTFparams.aberrations(:,3)*OTFparams.lambda; % change to length units
[OTF_orders,OTF_every_slide] = get_vectormodelOTF(OTFparams);              % compute vector model based OTF, SIM image sampling
% OTF_orders=abs(OTF_orders);
OTF_orders=gpuArray(single(real(OTF_orders)));
end
function [OTFinc,OTFinc2d_throughfocus] = get_vectormodelOTF(OTFparams)
%% Get pupil matrix
[~,~,wavevector,wavevectorzmed,~,PupilMatrix] = get_pupil_matrix(OTFparams);

%% Get field matrix
[XImage,YImage,ZImage,FieldMatrix] = get_field_matrix(PupilMatrix,wavevector,wavevectorzmed,OTFparams);

%% Get 3D PSF
PSF = get_psf(FieldMatrix,OTFparams);

%% Calculation of through-focus OTF for the focal stack by 2D-CZT in xy
[~,~,OTFinc2d_throughfocus] = get_throughfocusotf(PSF,XImage,YImage,OTFparams);

%% Calculation of 3D-OTF by 1D-CZT in z
[~,OTFinc] = get_otf3d(OTFinc2d_throughfocus,ZImage,OTFparams);

%% Masking out-of-band numerical noise to zero
OTFinc = do_OTFmasking3D(OTFinc,OTFparams);

%% Normalize the OTF
% centerpos = floor(size(OTFinc)/2)+1;
% if length(size(OTFinc))==3
%     OTFnorm = OTFinc(centerpos(1),centerpos(2),centerpos(3));
% end
% if length(size(OTFinc))==2
%     OTFnorm = OTFinc(centerpos(1),centerpos(2));
% end
OTFnorm=max(OTFinc(:));
OTFinc = OTFinc/OTFnorm;
end
function [XPupil,YPupil,wavevector,wavevectorzmed,Waberration,PupilMatrix] = get_pupil_matrix(OTFparams)

%% Basic parameter
NA = OTFparams.NA;
refmed = OTFparams.refmed;
refcov = OTFparams.refcov;
refimm = OTFparams.refimm;
refimmnom = OTFparams.refimmnom;
lambda = OTFparams.lambda;
Npupil = OTFparams.Npupil;

%% Feniel parameter
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = gpuArray(-PupilSize+DxyPupil/2:DxyPupil:PupilSize);
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
CosThetaMed = 1.0*sqrt(complex(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2));
CosThetaCov = 1.0*sqrt(complex(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2));
CosThetaImm = 1.0*sqrt(complex(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2));
CosThetaImmnom = 1.0*sqrt(complex(1-(XPupil.^2+YPupil.^2)*NA^2/refimmnom^2));
FresnelPmedcov = 2.0./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2.0./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2.0*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2.0*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP =  1.0*FresnelPmedcov.*FresnelPcovimm;
FresnelS =  1.0*FresnelSmedcov.*FresnelScovimm;

%% Entrance pupil function of plane
Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = 1.0*sqrt(complex(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2));
SinTheta = 1.0*sqrt(complex(1-CosTheta.^2));
pvec{1} = 1.0*FresnelP.*CosTheta.*CosPhi;
pvec{2} = 1.0*FresnelP.*CosTheta.*SinPhi;
pvec{3} = -1.0*FresnelP.*SinTheta;
svec{1} = -1.0*FresnelS.*SinPhi;
svec{2} = 1.0*FresnelS.*CosPhi;
svec{3} = 0;

%% XYZ direction component of P/S polarization
PolarizationVector = cell(2,3);
for jtel = 1:3
    PolarizationVector{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
    PolarizationVector{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
end

%% Aperture function
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
Amplitude = ApertureMask.*sqrt(CosThetaImm);

%% Correct phase difference
Waberration = 1.0*zeros(size(XPupil));
[zvals] = get_rimismatchpars(OTFparams);
Waberration = Waberration+zvals(1)*refimm*CosThetaImm-zvals(2)*refimmnom*CosThetaImmnom-zvals(3)*refmed*CosThetaMed;
PhaseFactor = exp(2*pi*1i*Waberration/lambda);

%% Entrance pupil matrix
PupilMatrix = cell(2,3);
for itel = 1:2
    for jtel = 1:3
        PupilMatrix{itel,jtel} = Amplitude.*PhaseFactor.*PolarizationVector{itel,jtel};
    end
end

%% Normalize
[normint_free] = get_normalization(PupilMatrix,OTFparams);
for itel = 1:2
    for jtel = 1:3
        PupilMatrix{itel,jtel} = PupilMatrix{itel,jtel}/sqrt(normint_free);
    end
end

%% Wave vector
wavevector = cell(1,3);
wavevector{1} = (2*pi*NA/lambda)*XPupil;
wavevector{2} = (2*pi*NA/lambda)*YPupil;
wavevector{3} = (2*pi*refimm/lambda)*CosThetaImm;
wavevectorzmed = (2*pi*refmed/lambda)*CosThetaMed;

end
function [zvals] = get_rimismatchpars(OTFparams)
% This function computes relevant parameters for the aberration function
% pertaining to refractive index mismatch. In particular, the routine
% computes the optimum z-stage position, given the nominal free working
% distance and imaging depth from the cover slip and the different
% refractive index values.
%
% relevant z-positions:
% zvals = [nominal stage position, free working distance, -image depth from cover slip]
%
% overall aberration function
% W = zvals_1*W_1 - zvals_2*W_2 - zvals_3*W_3
% W_1 = sqrt(refimm^2-rho^2*NA^2)
% W_2 = sqrt(refimmnom^2-rho^2*NA^2)
% W_3 = sqrt(refmed^2-rho^2*NA^2)
%
% The unknown parameter zvals_1 is computed by minimizing the variance in W
% across the pupil plane W_rms^2 = <W^2> - <W>^2, this is evaluated
% analyticaly with Mathematica, involving:
%
% Amat_{ij} = <W_i*W_j> - <W_i>*<W_j>
%
% minimizing W_rms^2 w.r.t. zvals_1 gives:
%
% zvals_1 = (Amat_{12}/Amat_{11})*zvals_2 + (Amat_{13}/Amat_{11})*zvals_3
%
% and a minimum/optimum rms wavefront aberration
%
% W_rms^2 = (Amat_{22} - Amat_{12}^2/Amat_{11})*zvals_2^2
%         + (Amat_{33} - Amat_{13}^2/Amat_{11})*zvals_3^2
%         + 2*(Amat_{23} - Amat_{12}*Amat_{13}/Amat_{11})*zvals_2*zvals_3
%
% The resulting expressions are numerically unstable in the paraxial regime,
% there a Taylor-series in NA is used. So far, it is assumed that all waves
% are propagating, i.e. NA < refmed in particular, so the outcome should be
% applied with caution to TIRF situations.
%
% copyright Sjoerd Stallinga, TU Delft, 2017

refins = [OTFparams.refimm OTFparams.refimmnom OTFparams.refmed];
zvals = [0 OTFparams.fwd -OTFparams.depth];
NA = OTFparams.NA;
% reduce NA in case of TIRF conditions
if (NA>OTFparams.refmed)
    NA = OTFparams.refmed;
end
paraxiallimit = 0.2;
K = length(refins);

if NA>paraxiallimit
    fsqav = zeros(K,1);
    fav = zeros(K,1);
    Amat = zeros(K,K);
    for jj = 1:K
        fsqav(jj) = refins(jj)^2-(1/2)*NA^2;
        fav(jj) = (2/3/NA^2)*(refins(jj)^3-(refins(jj)^2-NA^2)^(3/2));
        %     fav(jj) = (2/3)*(3*refins(jj)^4-3*refins(jj)^2*NA^2+NA^4)/(refins(jj)^3+(refins(jj)^2-NA^2)^(3/2));
        Amat(jj,jj) = fsqav(jj)-fav(jj)^2;
        for kk = 1:jj-1
            Amat(jj,kk) = (1/4/NA^2)*(refins(jj)*refins(kk)*(refins(jj)^2+refins(kk)^2)...
                -(refins(jj)^2+refins(kk)^2-2*NA^2)*sqrt(refins(jj)^2-NA^2)*sqrt(refins(kk)^2-NA^2)...
                +(refins(jj)^2-refins(kk)^2)^2*log((sqrt(refins(jj)^2-NA^2)+sqrt(refins(kk)^2-NA^2))/(refins(jj)+refins(kk))));
            Amat(jj,kk) = Amat(jj,kk)-fav(jj)*fav(kk);
            Amat(kk,jj) = Amat(jj,kk);
        end
    end
    zvalsratio = zeros(K,1);
    Wrmsratio = zeros(K,K);
    for jv = 2:K
        zvalsratio(jv) = Amat(1,jv)/Amat(1,1);
        for kv = 2:K
            Wrmsratio(jv,kv) = Amat(jv,kv)-Amat(1,jv)*Amat(1,kv)/Amat(1,1);
        end
    end
else
    % paraxial limit, Taylor-series in NA
    zvalsratio = zeros(K,1);
    Wrmsratio = [0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0];
    %     Wrmsratio = zeros(K,K);
    for jv = 2:K
        zvalsratio(jv) = refins(1)/refins(jv)+NA^2*(refins(1)^2-refins(jv)^2)/(4*refins(1)*refins(jv)^3);
        for kv = 2:K
            Wrmsratio(jv,kv) = NA^8*(refins(1)^2-refins(jv)^2)*(refins(1)^2-refins(kv)^2)/(11520*refins(1)^4*refins(jv)^3*refins(kv)^3);
        end
    end
end
zvals(1) = zvalsratio(2)*zvals(2)+zvalsratio(3)*zvals(3);


end
function [normint_free] = get_normalization(PupilMatrix,parameters)
% This function computes the PSF normalization by evaluating the energy
% flow through the lens aperture
%
% copyright Sjoerd Stallinga, TU Delft, 2017

% parameters
Npupil = parameters.Npupil;
NA = parameters.NA;
lambda = parameters.lambda;
pixelsize = parameters.pixelsize;

% intensity matrix
IntensityMatrix = zeros(3,3);
for itel = 1:3
    for jtel = 1:3
        for ztel = 1:2
            pupmat1 = PupilMatrix{ztel,itel};
            pupmat2 = PupilMatrix{ztel,jtel};
            IntensityMatrix(itel,jtel) = IntensityMatrix(itel,jtel)+...
                sum(sum(real(pupmat1.*conj(pupmat2))));
        end
    end
end

% normalization to take into account discretization correctly
% DxyPupil = 2*NA/lambda/Npupil;
% normfac = DxyPupil^2/(pixelsize)^2;
DxyPupil = 2/Npupil;
normfac = DxyPupil^2/(pixelsize*NA/lambda)^2;
IntensityMatrix = normfac*IntensityMatrix;

% evaluation normalization factors
normint_free = sum(diag(IntensityMatrix))/3;


end
function dataout = cztfunc(datain,A,B,D)
% This function evaluates the FT via the czt-algorithm
% arguments: datain = input data, dimensions K x N
%            A,B,D = auxiliary vectors computed in prechirpz, must have
%            lengths N, M, and L=N+M-1
% function value: dataout = output data, dimensions K x M

N = length(A);
M = length(B);
L = length(D);
K = size(datain,1);
Amt = repmat(A,K,1);
Bmt = repmat(B,K,1);
Dmt = repmat(D,K,1);

cztin =  zeros(K,L,'gpuArray');
cztin(:,1:N)= Amt.*datain;
temp = Dmt.*fft(cztin,[],2);
cztout = ifft(temp,[],2);
dataout = Bmt.*cztout(:,1:M);

end
function [XImage,YImage,ZImage,FieldMatrix] = get_field_matrix(PupilMatrix,wavevector,wavevectorzmed,OTFparams)

%% Basic parameter
NA = OTFparams.NA;
lambda = OTFparams.lambda;
xemit = OTFparams.xemit;
yemit = OTFparams.yemit;
zemit = OTFparams.zemit;
xrange = OTFparams.xrange;
yrange = OTFparams.yrange;
zmin = OTFparams.zrange(1);
zmax = OTFparams.zrange(2);
Npupil = OTFparams.Npupil;
Mx = OTFparams.Mx;
My = OTFparams.My;
Mz = OTFparams.Mz;
PupilSize = NA/lambda;
ImageSizex = xrange;
ImageSizey = yrange;
ImageSizez = (zmax-zmin)/2;
DxImage = 2*ImageSizex/Mx;
DyImage = 2*ImageSizey/My;
DzImage = 2*ImageSizez/Mz;
ximagelin = gpuArray(-ImageSizex+DxImage/2:DxImage:ImageSizex);
yimagelin =  gpuArray(-ImageSizey+DyImage/2:DyImage:ImageSizey);
[YImage,XImage] = meshgrid(yimagelin,ximagelin);
ZImage = gpuArray(zmin+DzImage/2:DzImage:zmax);

%% Aberration correction
Wpos = xemit*wavevector{1}+yemit*wavevector{2}+zemit*wavevectorzmed;

%% Auxiliary vector of DZT transformation
[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,Mx);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,My);

%% Field matrix
FieldMatrix = cell(2,3,Mz);
for jz = 1:numel(ZImage)
    zemitrun = ZImage(jz);
    PositionPhaseMask = exp(1i*(Wpos+zemitrun*wavevectorzmed));
    for itel = 1:2
        for jtel = 1:3
            PupilFunction = PositionPhaseMask.*PupilMatrix{itel,jtel};
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
        end
    end
end

end
function PSFout = do_pixel_blurring(PSFin,OTFparams)
% copyright Sjoerd Stallinga, TU Delft, 2018

oversampling = OTFparams.pixelsize/OTFparams.samplingdistance;
PSFsize = gpuArray(size(PSFin));
Mx = PSFsize(1);
My = PSFsize(2);
if mod(Mx,2)==1
    centerx = (Mx+1)/2;
else
    centerx = Mx/2+1;
end
if mod(My,2)==1
    centery = (My+1)/2;
else
    centery = My/2+1;
end

qxnorm = oversampling*((1:Mx)-centerx)/Mx;
qynorm = oversampling*((1:My)-centery)/My;
[Qx,Qy] = meshgrid(qxnorm,qynorm);
pixelblurkernel = sinc(Qx).*sinc(Qy);

if length(PSFsize)==2
    OTF = fftshift(fft2(PSFin));
    OTF = pixelblurkernel.*OTF;
    PSFout = ifft2(ifftshift(OTF));
    PSFout = real(PSFout);
end

if length(PSFsize)==3
    Mz = PSFsize(3);
    PSFout = zeros(Mx,My,Mz,'gpuArray');
    for jz = 1:Mz
        PSFslice = squeeze(PSFin(:,:,jz));
        OTF = fftshift(fft2(PSFslice));
        OTF = pixelblurkernel.*OTF;
        PSFout(:,:,jz) = ifft2(ifftshift(OTF));
        PSFout(:,:,jz) = real(PSFout(:,:,jz));
    end
end

end

function PSF = get_psf(FieldMatrix,OTFparams)

%% Average the field matrix to compute PSF
dims = size(FieldMatrix);
Mz = dims(3);
imdims = size(FieldMatrix{1,1,1});
Mx = imdims(1);
My = imdims(2);
PSF = zeros(Mx,My,Mz,'gpuArray');
for jz = 1:Mz
    for jtel = 1:3
        for itel = 1:2
            PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
        end
    end
end

%% Pixel blurring
PSF = do_pixel_blurring(PSF,OTFparams);

end
function [A,B,D] = prechirpz(xsize,qsize,N,M)
% This function evaluates the auxiliary vectors for the evaluation of the
% FT via the czt-algorithm
% arguments: xsize = window in real space abs(x)<xsize
%            qsize = window in Fourier space abs(q)<qsize
%            N = # sample points in real space (even)
%            M = # sample points in Fourier space (odd)
% function value: A,B,D = auxiliary vectors of lengths N, M, and L=N+M-1

L = N+M-1;
sigma = 2*pi*xsize*qsize/N/M;
Afac = exp(2*1i*sigma*(1-M));
Bfac = exp(2*1i*sigma*(1-N));
sqW = exp(2*1i*sigma);
W = sqW^2;
Gfac = (2*xsize/N)*exp(1i*sigma*(1-N)*(1-M));

Utmp = zeros(1,N,'gpuArray');
A = zeros(1,N,'gpuArray');
Utmp(1) = sqW*Afac;
A(1) = 1.0;
for i=2:N
    A(i) = Utmp(i-1)*A(i-1);
    Utmp(i) = Utmp(i-1)*W;
end

Utmp = zeros(1,M,'gpuArray');
B = ones(1,M,'gpuArray');
Utmp(1) = sqW*Bfac;
B(1) = Gfac;
for i=2:M
    B(i) = Utmp(i-1)*B(i-1);
    Utmp(i) = Utmp(i-1)*W;
end

Utmp = zeros(1,max(N,M)+1,'gpuArray');
Vtmp = zeros(1,max(N,M)+1,'gpuArray');
Utmp(1) = sqW;
Vtmp(1) = 1.0;
for i=2:max(N,M)+1
    Vtmp(i) = Utmp(i-1)*Vtmp(i-1);
    Utmp(i) = Utmp(i-1)*W;
end
D = ones(1,L,'gpuArray');
for i=1:M
    D(i) = conj(Vtmp(i));
end
for i=1:N
    D(L+1-i) = conj(Vtmp(i+1));
end

D = fft(D);

end
function [XSupport,YSupport,OTF] = get_throughfocusotf(PSF,XImage,YImage,OTFparams)

ImageSizex = (OTFparams.xrange);
ImageSizey = (OTFparams.yrange);
SupportSizex = (OTFparams.supportsizex);
SupportSizey = (OTFparams.supportsizey);
Nsupportx = (OTFparams.Nsupportx);
Nsupporty = (OTFparams.Nsupporty);
Mx = (size(PSF,1));
My = (size(PSF,2));
Mz = (size(PSF,3));

%% OTF support and sampling (in physical units)
DxSupport = 2*SupportSizex/Nsupportx;
DySupport = 2*SupportSizey/Nsupporty;
delqx = OTFparams.shiftsupport(1)*DxSupport;
delqy = OTFparams.shiftsupport(2)*DySupport;
xsupportlin = -SupportSizex+DxSupport/2:DxSupport:SupportSizex;
ysupportlin = -SupportSizey+DySupport/2:DySupport:SupportSizey;
[XSupport,YSupport] = meshgrid(xsupportlin-delqx,ysupportlin-delqy);

%% calculate auxiliary vectors for chirpz
[Ax,Bx,Dx] = prechirpz(ImageSizex,SupportSizex,Mx,Nsupportx);
[Ay,By,Dy] = prechirpz(ImageSizey,SupportSizex,My,Nsupporty);

%% calculation of through-focus OTF, and for each focus level the OTF peak is normalized to one
OTF = zeros(Nsupportx,Nsupporty,Mz,'gpuArray');
for jz = 1:Mz
    PSFslice = squeeze(PSF(:,:,jz));
    PSFslice = exp(-2*pi*1i*(delqx*XImage+delqy*YImage)).*PSFslice;
    IntermediateImage = transpose(cztfunc(PSFslice,Ay,By,Dy));
    tempim = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
    OTF(:,:,jz) = tempim/max(max(abs(tempim)));
end

end
function [ZSupport,OTF3D] = get_otf3d(OTFstack,ZImage,parameters)
% This function computes the 3D-OTF from a through-focus stack of the
% 2D-OTF by making an FT along the axial direction by applying a 1D-CZT
% along that direction.

%% Basic parameter
SupportSizez = parameters.supportsizez;
Nsupportz = parameters.Nsupportz;
zmin = parameters.zrange(1);
zmax = parameters.zrange(2);
ImageSizez = (zmax-zmin)/2;
[Nx,Ny,Mz] = size(OTFstack);

%% OTF support and sampling (in physical units)
DzSupport = 2*SupportSizez/Nsupportz;
delqz = parameters.shiftsupport(3)*DzSupport;
ZSupport = -SupportSizez+DzSupport/2:DzSupport:SupportSizez;
ZSupport = ZSupport-delqz;

%% calculate auxiliary vectors for chirpz
[A,B,D] = prechirpz(ImageSizez,SupportSizez,Mz,Nsupportz);

%% 1-D CZT
OTF3D = zeros(Nx,Ny,Nsupportz);OTFstack=gather(OTFstack,[]);ZImage=gather(ZImage);A=gather(A);B=gather(B);D=gather(D);
for ii = 1:Nx
    for jj = 1:Ny
        axialcut = squeeze(OTFstack(ii,jj,:))';
        axialcut = exp(-2*pi*1i*delqz*ZImage).*axialcut;
        OTF3D(ii,jj,:) = cztfunc1(axialcut,A,B,D);
    end
end
OTF3D(:,:,Nsupportz) = OTF3D(:,:,1);

%% Normalize
norm = max(abs(OTF3D(:)));
OTF3D = OTF3D/norm;

end
function OTFinc = do_OTFmasking3D(OTFinc,parameters)

%% Basic parameter
NA = parameters.NA;
lambda = parameters.lambda;
refmed = parameters.refmed;

%% OTF support and sampling (in physical units)
SupportSizex = parameters.supportsizex;
SupportSizey = parameters.supportsizey;
SupportSizez = parameters.supportsizez;
Nsupportx = parameters.Nsupportx;
Nsupporty = parameters.Nsupporty;
Nsupportz = parameters.Nsupportz;
Dqx = 2*SupportSizex/Nsupportx;
Dqy = 2*SupportSizex/Nsupporty;
Dqz = 2*SupportSizez/Nsupportz;
shiftqx = parameters.shiftsupport(1)*Dqx;
shiftqy = parameters.shiftsupport(2)*Dqy;
shiftqz = parameters.shiftsupport(3)*Dqz;
qx = -SupportSizex+Dqx/2:Dqx:SupportSizex;
qy = -SupportSizey+Dqy/2:Dqy:SupportSizey;
qz = -SupportSizez+Dqz/2:Dqz:SupportSizez;
qx = qx-shiftqx;
qy = qy-shiftqy;
qz = qz-shiftqz;
[spatfreqsx,spatfreqsy,spatfreqsz] = meshgrid(qx,qy,qz);

%% Masking
epsy = 1e2*eps;
q0 = refmed/lambda;
NAl = NA/lambda;
NBl = sqrt(q0^2-NAl^2);
spatfreqsxy = sqrt(spatfreqsx.^2+spatfreqsy.^2);
qcutoff = double(spatfreqsxy<=2*NAl).*(sqrt(complex(q0^2-(spatfreqsxy-NAl).^2))-NBl);
OTFmask = double(abs(spatfreqsz)<=qcutoff+Dqz/2);
OTFinc = OTFmask.*OTFinc;

end
function dataout = cztfunc1(datain,A,B,D)
% This function evaluates the FT via the czt-algorithm
% arguments: datain = input data, dimensions K x N
%            A,B,D = auxiliary vectors computed in prechirpz, must have
%            lengths N, M, and L=N+M-1
% function value: dataout = output data, dimensions K x M

N = length(A);
M = length(B);
L = length(D);
K = size(datain,1);
Amt = repmat(A,K,1);
Bmt = repmat(B,K,1);
Dmt = repmat(D,K,1);

cztin =  zeros(K,L);
cztin(:,1:N)= Amt.*datain;
temp = Dmt.*fft(cztin,[],2);
cztout = ifft(temp,[],2);
dataout = Bmt.*cztout(:,1:M);

end