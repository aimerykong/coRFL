% The code is used to compare receptive fields from different images (a
% pair as an example).
%
% vlfeat toolbox is borrowed here to support this code -- extracting SIFT
% features from images.
%
%
%
% Shu Kong (Aimery)
% aimerykong@gmail.com
% Feb. 2014

clear all;
close all;
clc;

addpath(genpath('./VLFEATROOT')); % vlfeat toolbox is required to be put here

%% parameters
sigma = .2; % sigma controls the Gaussian kernel that convert pair-wise distance to similarity
distanceRatio = 0.8; % distance ratio in comparing SIFT features by removing a part of features that are mainly negative ones

%% fetch two images
disp('Receptive Fields Compare for Image Pair...');

imNameA = 'superman3'; % 
imNameB = 'superman1'; % 

disp('SIFT extraction for image A...');
imgPath =  ['./data/', imNameA, '.jpg'];
[imFeaA, imA] = imSIFT(imgPath);
valjunk = sum(imFeaA.Fea.^2,1);
idx = find(valjunk > 0.7); % normalize the SIFT features
imFeaA.Fea(:,idx) = imFeaA.Fea(:,idx) ./ repmat( sqrt(valjunk(idx)), size(imFeaA.Fea,1), 1);

%[imFea] = SIFTinRF4Img_show(imgPath, flag_display);
disp('SIFT extraction for image B...');
imgPath =  ['./data/', imNameB, '.jpg'];
[imFeaB, imB] = imSIFT(imgPath);
valjunk = sum(imFeaB.Fea.^2,1);
idx = find(valjunk > 0.7); % normalize the SIFT features
imFeaB.Fea(:,idx) = imFeaB.Fea(:,idx) ./ repmat( sqrt(valjunk(idx)), size(imFeaB.Fea,1), 1);
clear idx valjunk

%% Receptive Fields generation
RFtemplates = genRFcandidates();
numRFperImg = size(RFtemplates, 2);

RFlist_A = RFtemplates; % adapt receptive field candidates according to the templates
RFlist_A = RFlist_A-1;
RFlist_A([1 3], :) = floor(RFlist_A([1 3], :)*(imFeaA.imsize(1)/100)) + 1;
RFlist_A([2 4], :) = floor(RFlist_A([2 4], :)*(imFeaA.imsize(2)/100)) + 1;

RFlist_B = RFtemplates; % adapt receptive field candidates according to the templates
RFlist_B = RFlist_B-1;
RFlist_B([1 3], :) = floor(RFlist_B([1 3], :)*(imFeaB.imsize(1)/100)) + 1;
RFlist_B([2 4], :) = floor(RFlist_B([2 4], :)*(imFeaB.imsize(2)/100)) + 1;

%% using a matrix to recods pairwise RF distance based on SIFT features
disp('pairwise distance of RFs in two images...');
aa = sum(imFeaA.Fea.^2, 1);
bb = sum(imFeaB.Fea.^2,1);
AB = repmat(aa', 1, size(imFeaB.Fea,2))+...
    repmat(bb, size(imFeaA.Fea,2), 1)-...
    2*imFeaA.Fea'*imFeaB.Fea;

SimilarityMat = zeros(numRFperImg, numRFperImg);
for i = 1:numRFperImg % from imageA
    idx = find( ...
        imFeaA.FeaLoc(1,:) >= RFlist_A(2,i) & ...
        imFeaA.FeaLoc(1,:) <= RFlist_A(4,i) & ...
        imFeaA.FeaLoc(2,:) >= RFlist_A(1,i) & ...
        imFeaA.FeaLoc(2,:) <= RFlist_A(3,i) ); % find the SIFT within the current RF
    RFsiftA = AB(idx,:);
    
    for j = 1:numRFperImg % from imageB
        idx = find( ...
            imFeaB.FeaLoc(1,:) >= RFlist_B(2,j) & ...
            imFeaB.FeaLoc(1,:) <= RFlist_B(4,j) & ...
            imFeaB.FeaLoc(2,:) >= RFlist_B(1,j) & ...
            imFeaB.FeaLoc(2,:) <= RFlist_B(3,j) ); 
        
        RFsiftB = RFsiftA(:,idx); % find the SIFT for comparison within the current RFs of imageA and imageB

        aa = min(RFsiftB,[],2);
        aa = aa(find(aa<=max(aa)*distanceRatio));
        bb = min(RFsiftB,[],1);
        bb = bb(find(bb<=max(bb)*distanceRatio));
        
        if isempty(aa) || isempty(bb)
            RFsiftB = 10; % no SIFT means this RF is invalide, thus allocating a large value to the pairwise distance
        else
            RFsiftB = sum(bb)/length(bb)+sum(aa)/length(aa);
        end
        SimilarityMat(i,j) = RFsiftB/2;
    end
end
clear aa bb AB; 
SimilarityMat = (SimilarityMat+SimilarityMat') / 2; % obtain the undirected graph (symmetric)

%%  display the similarity graph
SimilarityMat = exp(-SimilarityMat / sigma);
figure;
aa = max(SimilarityMat(:));
aa = find(SimilarityMat<aa*0.8);
SimilarityMat(aa) = 0;

imagesc(SimilarityMat);
title('similarity graph between RFs of two images');
ylabel('RF index of image-1');
xlabel('RF index of image-2');
