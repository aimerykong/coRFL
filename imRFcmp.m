function [SimilarityMat] = imRFcmp(pathIMGa, pathIMGb, RFtemplates, sigma, distanceRatio, parloc)
% The code is used to compare receptive fields from different images
%
% vlfeat toolbox is borrowed here to support this code -- extracting SIFT
% features from images.
%
%
%
% Shu Kong (Aimery)
% aimerykong@gmail.com
% Feb. 2014

%addpath(genpath('./VLFEATROOT')); % vlfeat toolbox is required to be put here

%% parameters
if nargin < 3
    RFtemplates = genRFcandidates();
end
if nargin < 4
    sigma = .5; % sigma controls the Gaussian kernel that convert pair-wise distance to similarity
end
if nargin < 5
    distanceRatio = .4; % distance ratio in comparing SIFT features by removing a part of features that are mainly negative ones
end
if nargin < 6
    parloc = 0.8; % penalize the locations of SIFT features
end

upperbound = 20;

%% fetch two images
imFeaA = imSIFT(pathIMGa);
imFeaA.sift = normalizeSIFT(imFeaA.sift); % normalize the SIFT features

%disp('SIFT extraction for image B...');
imFeaB = imSIFT(pathIMGb);
imFeaB.sift = normalizeSIFT(imFeaB.sift); % normalize the SIFT features

%% Receptive Fields generation
numRFperImg = size(RFtemplates, 2);

RFlist_A = RFtemplates; % adapt receptive field candidates according to the templates
RFlist_A([1 3], :) = floor( (RFlist_A([1 3], :)-1)*(imFeaA.imsize(1)/100))+1;
RFlist_A([2 4], :) = floor( (RFlist_A([2 4], :)-1)*(imFeaA.imsize(2)/100))+1;
% RFlist_A([1 3], :) = floor( (RFlist_A([1 3], :)*(imFeaA.imsize(1)/100)) );
% RFlist_A([2 4], :) = floor( (RFlist_A([2 4], :)*(imFeaA.imsize(2)/100)) );

RFlist_B = RFtemplates; % adapt receptive field candidates according to the templates
RFlist_B([1 3], :) = floor( (RFlist_B([1 3], :)-1)*(imFeaB.imsize(1)/100)) +1;
RFlist_B([2 4], :) = floor( (RFlist_B([2 4], :)-1)*(imFeaB.imsize(2)/100)) +1;
% RFlist_B([1 3], :) = floor( (RFlist_B([1 3], :)*(imFeaB.imsize(1)/100)) );
% RFlist_B([2 4], :) = floor( (RFlist_B([2 4], :)*(imFeaB.imsize(2)/100)) );

%% using a matrix to recods pairwise RF distance based on SIFT features
aa = sum(imFeaA.sift.^2, 1);
bb = sum(imFeaB.sift.^2,1);
AB = repmat(aa', 1, size(imFeaB.sift,2))+...
    repmat(bb, size(imFeaA.sift,2), 1)-...
    2*imFeaA.sift'*imFeaB.sift;

SimilarityMat = zeros(numRFperImg, numRFperImg);
SimilarityMat = checkRF(0.75, SimilarityMat, imFeaA.intraColor, imFeaB.intraColor, RFlist_A, RFlist_B);
colorSimMat = SimilarityMat;

for i = 1:numRFperImg % from imageA
    if sum(SimilarityMat(i,:)) == upperbound*size(SimilarityMat,2)
        continue;
    end
    idx = find( ...
        imFeaA.siftLoc(1,:) >= RFlist_A(2,i) & ...
        imFeaA.siftLoc(1,:) <= RFlist_A(4,i) & ...
        imFeaA.siftLoc(2,:) >= RFlist_A(1,i) & ...
        imFeaA.siftLoc(2,:) <= RFlist_A(3,i) ); % find the SIFT within the current RF
    RFsiftA = AB(idx,:);
    if isempty(idx)
        SimilarityMat(i,:) = upperbound; % a larger value to denote this pair is dissimilar
        continue;
    end
    
    RFsiftAloc = imFeaA.siftLoc(1:2, idx); % calculate SIFT location to achieve rough structures
    idx = [RFlist_A(2,i); RFlist_A(1,i)];
    RFsiftAloc = RFsiftAloc - repmat(idx, 1, size(RFsiftAloc,2));
    RFsiftAloc = RFsiftAloc ./ repmat([RFlist_A(4,i)-RFlist_A(2,i); RFlist_A(3,i)-RFlist_A(1,i)], 1, size(RFsiftAloc,2));
    RFsiftAloc = sqrt(sum(RFsiftAloc.^2,1));
    
   RFcandiA = imFeaA.Lab( RFlist_A(1,i):RFlist_A(3,i), RFlist_A(2,i):RFlist_A(4,i), 1:3 );
    %RFcandiA = imFeaA.im( RFlist_A(1,i):RFlist_A(3,i), RFlist_A(2,i):RFlist_A(4,i), 1:3 );
    
    for j = 1:numRFperImg % from imageB
        if sum(SimilarityMat(j,:)) == upperbound*size(SimilarityMat,2)
            continue;
        end
      
        idx = find( ...
            imFeaB.siftLoc(1,:) >= RFlist_B(2,j) & ...
            imFeaB.siftLoc(1,:) <= RFlist_B(4,j) & ...
            imFeaB.siftLoc(2,:) >= RFlist_B(1,j) & ...
            imFeaB.siftLoc(2,:) <= RFlist_B(3,j) ); 
        RFsiftB = RFsiftA(:,idx); % find the SIFT for comparison within the current RFs of imageA and imageB
        
        RFsiftBloc = imFeaB.siftLoc(1:2, idx); % calculate SIFT location to achieve rough structures
        idx = [RFlist_B(2,i); RFlist_B(1,i)];
        RFsiftBloc = RFsiftBloc - repmat(idx, 1, size(RFsiftBloc,2));
        RFsiftBloc = RFsiftBloc ./ repmat([RFlist_B(4,i)-RFlist_B(2,i); RFlist_B(3,i)-RFlist_B(1,i)], 1, size(RFsiftBloc,2));
        RFsiftBloc = sqrt(sum(RFsiftBloc.^2,1));
        
        RFcandiB = imFeaB.Lab( RFlist_B(1,i):RFlist_B(3,i), RFlist_B(2,i):RFlist_B(4,i), 1:3 );
        
        locGraph = repmat(RFsiftAloc', 1, length(RFsiftBloc)) + repmat(RFsiftBloc, length(RFsiftAloc), 1);
        
        RFsiftB = RFsiftB+parloc*locGraph;
        
        aa = min(RFsiftB,[],2);
        if isempty(aa)
            RFsiftB = upperbound;
        else
            a = find(aa<=max(aa)*distanceRatio);
            aa = aa(a);
            bb = min(RFsiftB,[],1);
            b = find(bb<=max(bb)*distanceRatio);
            bb = bb(b);
            if isempty(aa) || isempty(bb)
                RFsiftB = upperbound; % no SIFT means this RF is invalide, thus allocating a large value to the pairwise distance
            else
                RFsiftB = sum(bb)/length(bb)+sum(aa)/length(aa);              
            end
        end
        
        SimilarityMat(i,j) = RFsiftB/2;
        colorSimMat(i,j) = PEDcolor(RFcandiA, RFcandiB);%SimilarityMat(i,j);
    end
end
clear aa bb AB; 
SimilarityMat = (SimilarityMat+SimilarityMat') / 2; % obtain the undirected graph (symmetric)
colorSimMat = colorSimMat+colorSimMat';
colorSimMat = colorSimMat ./ max(colorSimMat(:));
SimilarityMat = SimilarityMat + 0.1*colorSimMat;

%%  display the similarity graph
SimilarityMat = exp(-SimilarityMat / sigma);
% [val, aa] = sort(SimilarityMat(:), 'descend');
% SimilarityMat( aa(round(length(aa)*0.3):end) ) = 0;
aa = max(SimilarityMat(:));
aa = find(SimilarityMat<aa*0.7);
SimilarityMat(aa) = 0;
