% The code is used to build similarity matrix for a set of images,
% recording their pairwise RF similarity.
%
% vlfeat toolbox is borrowed here to support this code -- extracting SIFT
% features from images.
%
% Current version does not adopt PED metric.
%
%
% Shu Kong (Aimery)
% aimerykong@gmail.com
% Feb. 2014

clear all;
close all;
clc;

addpath(genpath('./VLFEATROOT')); % vlfeat toolbox is required to be put here

%% data path and parameters (all data are from one common category)
CategName = 'cactus';
subPath = fullfile('./database4coRFL', CategName); %'./data/superman';%
RFtemplates = genRFcandidates('grid');
numRFperImg = size(RFtemplates, 2);

dirName = dir(subPath);
numImg = length(dirName)-2;

K = numImg*2; % select K desired receptive fields
tau = 2; %  requires tau>1 to control the pairwise similarity principle

%% build the graph by comparing pairwise RF distance of images
SimilarityMat = zeros(numImg*numRFperImg);

for a = 1:numImg
    pathIMGa = fullfile(subPath, dirName(a+2).name);
    for b = a+1:numImg
        disp(['build subgraph of RF distance for image pair: ' num2str(a) '&' num2str(b)]);
        pathIMGb = fullfile(subPath, dirName(b+2).name);
        [simMat] = imRFcmp(pathIMGa, pathIMGb, RFtemplates);

        SimilarityMat( 1+(a-1)*numRFperImg:a*numRFperImg, 1+(b-1)*numRFperImg:b*numRFperImg ) = simMat;
    end
end
SimilarityMat = (SimilarityMat+SimilarityMat')/2;

imagesc(SimilarityMat);
axis square;
title('similarity graph between RFs of images');
ylabel('RF index');
xlabel('RF index');

%% use the graph obtained above, feed it into the obj function to select the desirable RFs
lambda1 = 100; % balance term
lambda2 = .01; % center-bias
[RFpositive, imgIdx ] = CRFL4knowlege(SimilarityMat, tau, numImg, K, lambda1, lambda2, RFtemplates);
displayRF(RFpositive, imgIdx, subPath, RFtemplates);
