function [imgFea, im] = imSIFT(imgPath)
% given the path of one image, extract its SIFT features
%
%
%
%   Shu Kong (Aimery)
%   aimerykong@gmail.com
%   www.aimerykong.me
%   Feb. 2014

%% fetch image
im = imread( imgPath );

im = im2single( im ); % convert uint8 to double
im = rgb2gray(im);
imSize = size(im);
imSize = imSize(1:2);
if max(imSize) > 400
    im = imresize(im, 400/max(imSize)); % resize image with preserved aspect ratio
end
imgFea.imsize = size(im);

%% points of interest detection, and SIFT extraction
[imgFea.FeaLoc, imgFea.Fea] = vl_sift(im); 
imgFea.Fea = single(imgFea.Fea);
% A frame is a disk of center f(1:2), scale f(3) and orientation f(4)
% d: 128-dim descriptor


