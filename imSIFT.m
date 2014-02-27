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
if  min(size(imgPath)) == 1
    im = imread( imgPath );
else    
end

imSize = size(im);
imSize = imSize(1:2);
if max(imSize) > 300
    im = imresize(im, 300/max(imSize)); % resize image with preserved aspect ratio
end
im = im2single( im ); % convert uint8 to double
imgFea.im = im;
imgFea.imsize = size(im);

%% get the meaningful region w.r.t local contrast over color info
[L, a, b] = RGB2Lab( im );
imgFea.Lab(:,:,1) = L;
imgFea.Lab(:,:,2) = a;
imgFea.Lab(:,:,3) = b;

kernel = ones(24, 24);
sz = size(kernel);
step = floor(sz(1)/3);
kernel(  sz(1)/2-floor(step/2):sz(1)/2+floor(step/2), sz(2)/2-floor(step/2):sz(2)/2+floor(step/2) ) = -1;
idx = find(kernel==-1);
kernel(idx) = kernel(idx) / length(idx);
idx = find(kernel==1);
kernel(idx) = kernel(idx) / length(idx);

%L2 = abs(conv2(L, kernel, 'same'));
a2 = abs(conv2(a, kernel, 'same'));
b2 = abs(conv2(b, kernel, 'same'));
imgFea.intraColor = im2bw(a2+b2, .9);

%% points of interest detection, and SIFT extraction
im = rgb2gray(im);
[imgFea.siftLoc, imgFea.sift] = vl_sift(im, 'PeakThresh', 0.02, 'edgethresh', 30); 
imgFea.sift = single(imgFea.sift);

%{
close;
imshow(im);
h1   = vl_plotframe(imgFea.siftLoc) ; 
set(h1,'color','y','linewidth',2) ;
h2 = vl_plotsiftdescriptor(imgFea.sift, imgFea.siftLoc) ;
set(h2,'color','g','linewidth',1) ;
%}
