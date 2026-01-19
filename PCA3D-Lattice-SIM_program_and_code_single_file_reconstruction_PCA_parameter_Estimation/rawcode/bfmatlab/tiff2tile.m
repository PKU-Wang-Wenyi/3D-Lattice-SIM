% Tiling gray-scale tif files into a single large gray-scale tif file,
% using a tiled color image made by Image Composite Editor (ICE) as a
% reference.
% 
% Put following image files into a single folder:
%   - Tiled color image made by ICE (_stitch.tiff)
%   - Non-tiled gray-scale images (B/G/R.tif)
% 
% Takuma Kitanishi, OCU, 2011/11/02

clear;
fclose all;
figure(1); clf

% ----- indicate data folder -----
basepath = 'D:\data\from\Umaba\171101\171024\tile'

% get file list
fileT = dir(fullfile(basepath,'*_stitch.tiff'));
fileI = dir(fullfile(basepath,'*I.tif'));
fileR = dir(fullfile(basepath,'*R.tif'));
fileG = dir(fullfile(basepath,'*G.tif'));
fileB = dir(fullfile(basepath,'*B.tif'));

% input check
if length(fileT)~=1
    error('ERROR: Put a tiled image (_stitch.tiff) in the folder!')
end
if ~isequal(length(fileI),length(fileR),length(fileG),length(fileB))
    error('ERROR: Numbers of gray-scale images are unequal!')
end

% tiled color image
ImT = imread(fullfile(basepath,fileT.name));

% original tiled image, blue channel
ImTB = ImT(:,:,3);
H1 = size(ImTB,1);
W1 = size(ImTB,2);
% 
subplot(2,3,[1 2])
imagesc(ImTB); hold on
colormap default
axis image
title(fileT.name)

% resized tiled image, blue channel
scale = 0.1;
imTB = imresize(ImTB,scale);
h1 = size(imTB,1);
w1 = size(imTB,2);
siz2 = 1024;     % image size for 2nd alignment

IJ = nan(length(fileB),2);
for ii=1:length(fileB)
    
    % original gray image
    ImB = imread(fullfile(basepath,fileB(ii).name));
    H2 = size(ImB,1);
    W2 = size(ImB,2);
    %
    subplot(2,3,3)
    imagesc(ImB)
    title(fileB(ii).name)
    axis image
    
    % resized gray image
    imB = imresize(ImB,scale);
    h2 = size(imB,1);
    w2 = size(imB,2);
    
    % 1st alignment
    r = nan(h1-h2+1,w1-w2+1);
    for m=1:h1-h2+1
        for n=1:w1-w2+1
%             r(m,n) = corr2(imTB(m:m+h2-1,n:n+w2-1),imB);
            a = imTB(m:m+h2-1,n:n+w2-1);
            b = imB;
            a = single(a(:)-mean(a(:)));
            b = single(b(:)-mean(b(:)));
            r(m,n) = a'*b/(norm(a)*norm(b));
        end
    end
    [~,ind] = max(r(:));
    [i,j] = ind2sub(size(r),ind);
    i = round(i/scale);
    j = round(j/scale);
    %
    subplot(2,3,[4 5])
    imagesc(r)
    caxis([0 1])
    colorbar
    axis image
    title('1st alignment')
    
    % 2nd alignment
    R = nan(20);
    for m=1:20
        for n=1:20
            if i+m-20>0 && i+m+siz2-21<=H1 && j+n-20>0 && j+n+siz2-21<=W1
                a = ImTB(i+m-20:i+m+siz2-21,j+n-20:j+n+siz2-21);
                b = ImB(1:siz2,1:siz2);
                a = single(a(:)-mean(a(:)));
                b = single(b(:)-mean(b(:)));
                R(m,n) = a'*b/(norm(a)*norm(b));
            end
        end
    end
    [~,IND] = max(R(:));
    [I,J] = ind2sub(size(R),IND);
    %
    subplot(2,3,6)
    imagesc(R)
    axis image
    colorbar
    title('2nd alignment')
    
    % alignmennt coordinate
    IJ(ii,:) = [i+I-20 j+J-20];
    % 
    figure(1); 
    subplot(2,3,[1 2])
    rectangle('Position',[IJ(ii,2) IJ(ii,1) W2-1 H2-1 ],'EdgeColor','r')
    drawnow;
end

% construct images
tileB = uint8(zeros(size(ImTB)));
tileG = uint8(zeros(size(ImTB)));
tileR = uint8(zeros(size(ImTB)));
tileI = uint8(zeros(size(ImTB)));
for ii=1:length(fileB)
    tileB(IJ(ii,1):IJ(ii,1)+H2-1,IJ(ii,2):IJ(ii,2)+W2-1) = imread(fullfile(basepath,fileB(ii).name));
    tileG(IJ(ii,1):IJ(ii,1)+H2-1,IJ(ii,2):IJ(ii,2)+W2-1) = imread(fullfile(basepath,fileG(ii).name));
    tileR(IJ(ii,1):IJ(ii,1)+H2-1,IJ(ii,2):IJ(ii,2)+W2-1) = imread(fullfile(basepath,fileR(ii).name));
    tileI(IJ(ii,1):IJ(ii,1)+H2-1,IJ(ii,2):IJ(ii,2)+W2-1) = imread(fullfile(basepath,fileI(ii).name));
end

% save images
imwrite(tileB,[basepath '\' fileT.name(1:end-5) 'B.tiff'])
imwrite(tileG,[basepath '\' fileT.name(1:end-5) 'G.tiff'])
imwrite(tileR,[basepath '\' fileT.name(1:end-5) 'R.tiff'])
imwrite(tileI,[basepath '\' fileT.name(1:end-5) 'I.tiff'])


    