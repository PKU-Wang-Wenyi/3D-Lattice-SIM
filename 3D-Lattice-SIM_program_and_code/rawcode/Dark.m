function result_final = Dark(src,backgroundremoval)
NA=1.49;
emwavelength=610;
pixelsize = 65;
% backgroundremoval=0.175;
% dst,src,backgroundremoval,NA,emwavelength,pixelsize
%% Read data
image0 = gpuArray((double((src))));
% image0=importImages(image0);
image0 = 255*(image0 - min(min(image0)))./(max(max(image0))-min(min(image0)));
[Nx0,Ny0,~] = size(image0);
temp=image0;
temp(:,1,:)=image0(:,2,:);
temp(:,end,:)=image0(:,end-1,:);
image0=temp;
% figure;imshow(image0(:,:,1),[])
[Nx,Ny,~] = size(image0);
if Ny>Nx
    image0(Nx+1:Ny,:,:)=0;
elseif Ny<Nx
    image0(:,Ny+1:Nx,:)=0;
end
[Nx,Ny,Nz] = size(image0);

%% Reconstruction parameters
backgroundremoval=backgroundremoval*backgroundremoval;
if backgroundremoval<0.6
    background = 0; % 0-middle,1-severve
else
    background = 1;
end
pad = 1;        %1-sysemtic,0-pad0
denoise = 0;    % Guassion denoise
thres = 150;     % Threshold to distinguish background and information
divide = 0.5;

%% Padding edge
pad_size = 100;
result_stack = gpuArray.zeros(Nx,Ny,Nz,'double');
Lo_process_stack = gpuArray.zeros(Nx,Ny,Nz,'double');
Hi_stack = gpuArray.zeros(Nx,Ny,Nz,'double');
image=gpuArray.zeros(Nx+2*(floor(Nx/pad_size)+1),Ny+2*(floor(Ny/pad_size)+1),Nz,'double');
for jz = 1:Nz
    if pad ==1
        image(:,:,jz) = padarray(image0(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1],'symmetric');
    else
        image(:,:,jz) = padarray(image0(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1]);
    end
end

%% Basic parameters
[params.Nx,params.Ny,~] = size(image);
params.NA = NA;
params.emwavelength = emwavelength;
params.pixelsize = pixelsize;
params.factor = backgroundremoval;

%% Background setting
if background == 1
    maxtime = 2;
    deg_matrix = [6,4,1.2];   % 3-10
    dep_matrix = [3,3,2];   % 0.7-2
    hl_matrix = [1,1,1];    % 3-8
elseif background == 0
    maxtime=1;
    deg_matrix = 6;   % 3-10
    dep_matrix = 2;   % 0.7-2
    hl_matrix = 1;    % 3-8
end

%% Dark sectioning
for time = 1:maxtime
    for jz = 1:Nz
        deg = deg_matrix(maxtime);   % 3-10
        dep = dep_matrix(maxtime);   % 0.7-2
        hl = hl_matrix(maxtime);    % 3-8
        % Seperate spectrum and confirm block size
        [Hi,Lo,lp,EL] = separateHiLo(image(:,:,jz),params,deg,divide);
        block_size = confirm_block(params,lp);
        % Remove background for low-frequency part
        Lo_process = dehaze_fast2(Lo, 0.95, block_size, EL,dep,thres);
        result = Lo_process/hl + Hi;
        % Cutting edge
        Lo_process = Lo_process(floor(Nx/pad_size)+2:floor(Nx/pad_size)+Nx+1,floor(Ny/pad_size)+2:floor(Ny/pad_size)+Ny+1);
        % Lo = Lo(floor(Nx/pad_size)+2:floor(Nx/pad_size)+Nx+1,floor(Ny/pad_size)+2:floor(Ny/pad_size)+Ny+1);
        Hi = Hi(floor(Nx/pad_size)+2:floor(Nx/pad_size)+Nx+1,floor(Ny/pad_size)+2:floor(Ny/pad_size)+Ny+1);
        result = result(floor(Nx/pad_size)+2:floor(Nx/pad_size)+Nx+1,floor(Ny/pad_size)+2:floor(Ny/pad_size)+Ny+1);
        % Saving results
        result_stack(:,:,jz) = result;
        Lo_process_stack(:,:,jz) = Lo_process;
        Hi_stack(:,:,jz) = Hi;
    end
    image0 = result_stack;
    for jz = 1:Nz
        if pad==1
            image(:,:,jz) = padarray(image0(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1],'symmetric');
        else
            image(:,:,jz) = padarray(image0(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1]);
        end
    end
end

%% Denoise
result_final = zeros(Nx,Ny,Nz);
for jz = 1:Nz
    if pad ==1
        temp = padarray(result_stack(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1],'symmetric');
    else
        temp = padarray(result_stack(:,:,jz),[floor(Nx/pad_size)+1,floor(Ny/pad_size)+1]);
    end
    if denoise == 0
        temp1 = temp;
    else
        W = fspecial('gaussian',[2,2],1);
        temp1 = imfilter(temp, W, 'replicate');
    end
    result_final(:,:,jz) = temp1(floor(Nx/pad_size)+2:floor(Nx/pad_size)+Nx+1,floor(Ny/pad_size)+2:floor(Ny/pad_size)+Ny+1);
end

%% Select regin of interst

if Nx0~=Nx || Ny0~=Ny
    % if Nx>Nx0
    %     result_final(Nx0+1:Nx,:,:)=[];
    % end
    % if Ny0>Ny
    %     result_final(:,Ny0+1:Ny,:)=[];
    % end
    temp=result_final;
    clear result_final
    result_final=temp(1:Nx0,1:Ny0,:);
end

%% Saving results
% % maxnum = max(max(max(result_final)));
% final_image = uint16(result_final);
% stackfilename = dst;
% if exist(dst,'file')
%     delete(dst)
% end
% 
% for jj=1:size(final_image,3)
%     imwrite(final_image(:,:,jj),dst,'writemode','append')
% end
result_final=uint16((result_final));
% stackfilename = dst;
% imstackwrite(result_final, stackfilename)
end
function imstackwrite(stack, filename)
% Open the TIFF file in 'w8' mode to create a BigTIFF file
im = Tiff(filename, 'w8');

% Update infostruct for 16-bit data
infostruct.ImageLength = size(stack, 1);
infostruct.ImageWidth = size(stack, 2);
infostruct.Photometric = Tiff.Photometric.MinIsBlack;
infostruct.BitsPerSample = 16; % Change to 16 bits per sample
infostruct.SampleFormat = Tiff.SampleFormat.UInt; % Keep as unsigned integer
infostruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

for k = 1:size(stack, 3)
    im.setTag(infostruct);
    im.write(uint16(stack(:, :, k))); % Write as 16-bit unsigned integer
    im.writeDirectory();
end

im.close();

end

function sum_img = window_sum_filter(image, r)

% sum_img(x, y) = = sum(sum(image(x-r:x+r, y-r:y+r)));

[h, w] = size(image);
% sum_img = zeros(size(image));
sum_img = gpuArray.zeros(size(image));
% Y axis
im_cum = cumsum(image, 1);

sum_img(1:r+1, :) = im_cum(1+r:2*r+1, :);
sum_img(r+2:h-r, :) = im_cum(2*r+2:h, :) - im_cum(1:h-2*r-1, :);
sum_img(h-r+1:h, :) = repmat(im_cum(h, :), [r, 1]) - im_cum(h-2*r:h-r-1, :);

% X axis
im_cum = cumsum(sum_img, 2);

sum_img(:, 1:r+1) = im_cum(:, 1+r:2*r+1);
sum_img(:, r+2:w-r) = im_cum(:, 2*r+2:w) - im_cum(:, 1:w-2*r-1);
sum_img(:, w-r+1:w) = repmat(im_cum(:, w), [1, r]) - im_cum(:, w-2*r:w-r-1);

end
function [Hi,Lo,lp,EL] = separateHiLo(image,params,deg,divide)

%% 基本参数
Nx = params.Nx;
Ny = params.Ny;
NA = params.NA;
emwavelength = params.emwavelength;
pixel_size = params.pixelsize;

%% 其他参数
% res = 0.5 * emwavelength / NA/ params.factor;     % resolution
res = 0.5 * emwavelength / NA;
k_m = Ny / (res / pixel_size(1));    % objective cut-off frequency ???
kc = nearest(k_m * 0.2);             % cut-off frequency between hp and lp filter
sigmaLP = kc*2/2.355;                % Finding sigma value for low pass

%% 滤波
lp = lpgauss(Nx,Ny,sigmaLP*2*divide);
hp = hpgauss(Nx,Ny,sigmaLP*2*divide);
elp = lpgauss(Nx,Ny,sigmaLP/deg);
% ehp = hpgauss(Nx,Ny,sigmaLP/deg);
% figure;mesh(lp)
%% 得到高低频和极低频率
% Hi = real(ifft2(fft2(image).*hp));
% Lo = real(ifft2(fft2(image).*lp));

fft_image = fftshift(fft2(image));
Hi = real(ifft2(ifftshift(fft_image.*fftshift(hp))));
Lo = real(ifft2(ifftshift(fft_image.*fftshift(lp))));

% elp = zeros(Nx,Ny);
% elp(floor((Nx+1)/2)-1:floor((Nx+1)/2)+1,floor((Nx+1)/2)-1:floor((Nx+1)/2)+1) = 1;
EL = real(ifft2(fft2(image).*elp));
% EH = real(ifft2(fft2(image).*ehp));
end
function PSF = PSF_Generator(lambada,pixelsize,NA,w,factor)
w=(double(w));
[X,Y]=meshgrid(linspace(0,w-1,w),linspace(0,w-1,w));
scale=2*pi*NA/lambada*pixelsize;
scale=scale*factor;

R=sqrt(min(X,abs(X-w)).^2+min(Y,abs(Y-w)).^2);
PSF=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2;
PSF = PSF/(sum(sum(PSF)));
PSF=fftshift(PSF);

end
function [ out ] = hpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.
out=1-lpgauss(H,W,SIGMA);
end
function [ out ] = lpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.
H = gpuArray(double(H));
W = gpuArray(double(W));
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
% temp0 = -floor(W/2);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
temp = -(x.^2/(kcx^2)+y.^2/(kcy^2));
out = ifftshift(exp(temp));
% out = ifftshift(exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));
end
function ImStack = imstackread(img)
% This function can import 3D image stack.
% The bit depth is limited to 8,16,32-bit
% Written by Ethan Zhao, Jan. 2021
% Tutorial: https://zhuanlan.zhihu.com/p/326688179/

imgInfo = imfinfo(img);
imgRow = imgInfo(1).Height;
imgCol = imgInfo(1).Width;
imgDepth = length(imgInfo);
imgBitDepth = ['uint',num2str(imgInfo(1).BitDepth)];
ImStack = zeros(imgRow, imgCol, imgDepth,imgBitDepth);
for ii = 1 : length(imgInfo)
    ImStack(:,:,ii) = imread(img, ii);
end

end
function q = guided_filter(guide, target, radius, eps)
% Copyright (c) 2014 Stephen Tierney

[h, w] = size(guide);

avg_denom = window_sum_filter(gpuArray.ones(h, w,'double'), radius);

mean_g = window_sum_filter(guide, radius) ./ avg_denom;
mean_t = window_sum_filter(target, radius) ./ avg_denom;

corr_gg = window_sum_filter(guide .* guide, radius) ./ avg_denom;
corr_gt = window_sum_filter(guide .* target, radius) ./ avg_denom;

var_g = corr_gg - mean_g .* mean_g;
cov_gt = corr_gt - mean_g .* mean_t;

a = cov_gt ./ (var_g + eps);
b = mean_t - a .* mean_g;

mean_a = window_sum_filter(a, radius) ./ avg_denom;
mean_b = window_sum_filter(b, radius) ./ avg_denom;

q = mean_a .* guide + mean_b;

end
function trans_est = get_transmission_estimate(rep_atmosphere, image, omega, win_size)
% Copyright (c) 2014 Stephen Tierney
% [m, n, ~] = size(image);

%rep_atmosphere = repmat(atmosphere, m, n);

trans_est = 1 - omega * get_dark_channel( image ./ rep_atmosphere, win_size);

end
function radiance = get_radiance(rep_atmosphere,image, transmission)
% Copyright (c) 2014 Stephen Tierney

% [m, n, ~] = size(image);

max_transmission = max(transmission, 0.1);

radiance = ((image - rep_atmosphere) ./ max_transmission) + rep_atmosphere;

end



function dark_channel = get_dark_channel(image, win_size)
% [m, n, ~] = size(image);
% pad_size = floor(win_size/2);
% padded_image = padarray(image, [pad_size pad_size], Inf);
% dark_channel = zeros(m, n);
% for j = 1 : m
%     for i = 1 : n
%         patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
%         dark_channel1(j,i) = min(patch(:));
%     end
% end
dark_channel = movmin((image), win_size);
end
function atmosphere = get_atmosphere(image, dark_channel)
% Copyright (c) 2014 Stephen Tierney
% image=gather(image);
% dark_channel=gather(dark_channel);
[m, n, ~] = size(image);
n_pixels = m * n;

n_search_pixels = floor(n_pixels * 0.01);

dark_vec = reshape(dark_channel, n_pixels, 1);

image_vec = reshape(image, n_pixels,1);

[~, indices] = sort(dark_vec, 'descend');

% accumulator = 0;
%
% for k = 1 : n_search_pixels
%     accumulator = accumulator + image_vec(indices(k),:);
% end
accumulator=sum(image_vec(indices(1 : n_search_pixels),:),'all');
atmosphere =(accumulator / n_search_pixels) ;

end
function [ radiance ] = dehaze_fast2(  image, omega, win_size,EL,dep,thres )
% Copyright (c) 2014 Stephen Tierney


[Nx,Ny] = size(image);
if ~exist('omega', 'var')
    omega = 0.95;
end

if ~exist('win_size', 'var')
    win_size = 15;
end

r = 15;
res = 0.001;

[m, n, ~] = size(image);

Mask = zeros(Nx,Ny,'double');
Mask(image<thres) = 1;
dark_channel = get_dark_channel(image.*Mask, win_size);
min_atmosphere = get_atmosphere(image.*Mask, dark_channel);

dark_channel = get_dark_channel(image, win_size);
max_atmosphere = get_atmosphere(image, dark_channel);
EL = EL - min(min(EL));
rep_atmosphere_process = EL/max(max(EL))*(max_atmosphere-min_atmosphere)+min_atmosphere;
rep_atmosphere_process = dep*rep_atmosphere_process;
trans_est = get_transmission_estimate(rep_atmosphere_process,image, omega, win_size);
x = guided_filter(image, trans_est, r, res);
transmission = reshape(x, m, n);
radiance = get_radiance(rep_atmosphere_process,image, transmission);


end
function block_size = confirm_block(params,lp)
PSF = PSF_Generator(params.emwavelength,params.pixelsize,params.NA,params.Nx,params.factor);
PSF_Lo = abs(ifft2(fftshift(fft2(PSF)).*fftshift(lp)));
PSF_Lo = PSF_Lo./max(max(PSF_Lo));
% figure;plot(PSF(floor(params.Nx/2)+1,:))
for count_x = floor(params.Nx/2):params.Nx
    if PSF_Lo(count_x,floor(params.Nx/2)) <0.01
        break;
    end
end
block_size = count_x-floor(params.Nx/2);
end

