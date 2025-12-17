% Batch convert from czi files to color-separated tiff files.
% 
% Requires Bio-formats toolbox.
% https://www.openmicroscopy.org/bio-formats/
% 
% Takuma Kitanishi, OCU, 2017/11/01

clear all;
close all;
fclose all

% ----- indicate the folder containing czi files -----
basepath = 'D:\data\from\Umaba\171101\171024\czi'

% get czi file list
filelist = dir(fullfile(basepath,'*.czi'));

for ii=1:length(filelist)
    % open a czi file
    filename = filelist(ii).name;
    im = bfopen(fullfile(basepath,filename));
    
    % split color
    I = im{1,1}{1,1};
    R = im{1,1}{2,1};
    G = im{1,1}{3,1};
    B = im{1,1}{4,1};
    
    % show image
    figure(1)
    subplot(221); imagesc(I); title('IR');
    subplot(222); imagesc(R); title('R');
    subplot(223); imagesc(G); title('G');
    subplot(224); imagesc(B); title('B');
    drawnow;
    
    % save
    imwrite(I,[basepath '\' filename(1:end-4) 'I.tif'])
    imwrite(R,[basepath '\' filename(1:end-4) 'R.tif'])
    imwrite(G,[basepath '\' filename(1:end-4) 'G.tif'])
    imwrite(B,[basepath '\' filename(1:end-4) 'B.tif'])
end
