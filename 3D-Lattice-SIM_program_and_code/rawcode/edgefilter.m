function [out]=edgefilter(images)
if ndims(images)==4
    N = size(images, 3)*size(images, 4);
else
    N = size(images, 3)*size(images, 4)*size(images, 4);
end
L = size(images, 1);
if L > 256
    px = 10;
else
    px = 0;
end

for I = 1:N
    images(:,:,I) = fadeBorderCos(images(:,:,I), px);
end
out = images;
end
function [out] = fadeBorderCos(img, px)
[h, w] = size(img);
fac = 1 / px * pi / 2;
% 只有当px>0时才处理边缘淡化效果
if px > 0
    % 计算淡化权重
    weight = sin((0:px-1) * fac).^2;
    weightTopBottom = repmat(weight, w, 1).';
    weightLeftRight = repmat(weight, h, 1);
    % 应用权重
    img(1:px, :) = img(1:px, :) .* weightTopBottom;
    img(end-px+1:end, :) = img(end-px+1:end, :) .* flipud(weightTopBottom);
    img(:, 1:px) = img(:, 1:px) .* weightLeftRight;
    img(:, end-px+1:end) = img(:, end-px+1:end) .* fliplr(weightLeftRight);
end
out = img;
end