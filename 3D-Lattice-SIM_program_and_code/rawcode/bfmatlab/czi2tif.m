% Batch convert from czi files to color-merged tiff files.
% 
% Requires Bio-formats toolbox.
% https://www.openmicroscopy.org/bio-formats/
% 
% Takuma Kitanishi, OCU, 2017/11/01

clear all;
close all;
fclose all

% ----- indicate the folder containing czi files -----
basepath = 'D:\data\imLSM700\171102-SUBantero\tk0057\AVITN'
cell2ch(1) = 2;
cell2ch(2) = 3;

% get czi file list
filelist = dir(fullfile(basepath,'*.czi'));

for ii=1:length(filelist)
    % open a czi file
    filename = filelist(ii).name;
    im = bfopen(fullfile(basepath,filename));
    
    % number of channels
    n = length(im{1,1});
    
    % construct an image
    tif = uint8(zeros([size(im{1,1}{1,1}) 3]));
    for jj=1:n
        tif(:,:,cell2ch(jj)) = uint8((im{1,1}{jj,1}+1)/16-1);
    end
    
    % show image
    figure(1)
    image(tif)
    axis image
    title(filename)
    drawnow;
    
    % save
    imwrite(tif,[basepath '\' filename(1:end-4) '.tif'])
end
