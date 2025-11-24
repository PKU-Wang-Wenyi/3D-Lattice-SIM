function imwritestack_16(stack, filename)
stack=uint16(stack);
var_info = whos('stack');
size_in_GB = var_info.bytes / (1000^3);
if size_in_GB>3.99
    t = Tiff(filename, 'w8');
else
    t = Tiff(filename, 'w');
end
tagstruct.ImageLength = size(stack, 1);
tagstruct.ImageWidth = size(stack, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
for k = 1:size(stack, 3)
    t.setTag(tagstruct)
    t.write((stack(:, :, k)));
    t.writeDirectory();
end
t.close();
end