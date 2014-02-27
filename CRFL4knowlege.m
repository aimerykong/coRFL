function [RFpositive, imgIdx ]= CRFL4knowlege(SimilarityGraph, tau, numImg, K, lambda1, lambda2, RFtemplates)
%
%
%   Shu Kong (Aimery)
%   aimerykong@gmail.com
%   www.aimerykong.me
%   Dec. 2013


%% prelude
if nargin < 7
    RFtemplates = genRFcandidates();
end
numRFperImg = size(RFtemplates, 2);

imgIdx = zeros(1, numImg*numRFperImg);
for count = 1:numImg
    imgIdx( 1+numRFperImg*(count-1):count*numRFperImg ) = count;
end

%% create the center-bias constraint (to maximize)
CenterBiasFavor = RFtemplates(end,:);
CenterBiasFavor = CenterBiasFavor + mean(CenterBiasFavor);
CenterBiasFavor = sqrt(CenterBiasFavor/10);
CenterBiasFavor = CenterBiasFavor + 1;
CenterBiasFavor = 1 ./ CenterBiasFavor;
CenterBiasFavor = CenterBiasFavor ./ max(CenterBiasFavor);
CenterBiasFavor = repmat( CenterBiasFavor, 1, numImg );

%% greedily optimize the submodular function and select the most desired RF's
fprintf('\tseek for %d most desirable RF candidates... ', K);
SelectedFlag = zeros( 1, size(SimilarityGraph, 2) );

C = 1+tau*sum(SimilarityGraph(:));

CountNumRFperImg = zeros(1, numImg);
RFpositive = zeros(1, K);

flag = true;
k = 1;
while flag
    fprintf('.');
    
    SelectAlready = find(SelectedFlag == 1);
    Hcur = C + sum(sum( SimilarityGraph(SelectAlready, :) ));
    SelectNotYet = find(SelectedFlag ~= 1);
    Hcur = Hcur - tau*sum(sum( SimilarityGraph(SelectNotYet, :) ));
    
    Hcur = log(Hcur) + lambda1*sum(log(1+CountNumRFperImg)) + lambda2*sum(CenterBiasFavor(SelectAlready));
    %Hcur = log(Hcur) + lambda1*sum(log(1+CountNumRFperImg));
    
    predictH = zeros(1, length(SelectNotYet));
    g_balance = repmat( CountNumRFperImg(:), 1, size(predictH,2) );
    a = ( 0:length(SelectNotYet)-1 ) * size( g_balance, 1 );
    g_balance(a+imgIdx(SelectNotYet)) = g_balance(a+imgIdx(SelectNotYet)) + 1;
    
    
    PredictSelectSet = [repmat( SelectAlready(:), 1, length(SelectNotYet) ); SelectNotYet];
    PredictExclusiveSet = repmat(SelectNotYet(:), 1, length(SelectNotYet(:)) );
    PredictExclusiveSet = PredictExclusiveSet - diag( diag(PredictExclusiveSet) );
    PredictExclusiveSet = PredictExclusiveSet + diag( PredictExclusiveSet(1,:) );
    PredictExclusiveSet = PredictExclusiveSet(2:end,:);
    
    H_A = repmat( C + sum(sum(SimilarityGraph( SelectAlready, : ))),  1, length(SelectNotYet) );
    
    tmpGain = sum(SimilarityGraph(SelectNotYet(:), :), 2);
    H_A = H_A(:) + tmpGain;
    minusH_A = sum(sum(   SimilarityGraph( SelectNotYet, : )    ));
    minusH_A = minusH_A - tmpGain;
    H_A = H_A - tau*minusH_A;
    
    predictH = log( H_A' ) + lambda1*sum( log(1+g_balance), 1 ) ...
        + lambda2*( sum(CenterBiasFavor(SelectAlready))+ CenterBiasFavor(SelectNotYet) );
    
    %% calculate the marginal gain and find the most desirable exemplar
    predictGain = predictH(:) - Hcur; % calculate info gain by openning a specific facility
    [y, idx] = max(predictGain); % find the most informative facility which brings the most gain
    
    if y <= 0 || k >= K % jump out loop when no gains or fixed iterations exceeded
        SelectedFlag(SelectNotYet(idx)) = 1; % select the most informative one
        CountNumRFperImg(imgIdx(SelectNotYet(idx))) = CountNumRFperImg(imgIdx(SelectNotYet(idx))) + 1;
        RFpositive(k) = SelectNotYet(idx);
        flag = false;
        break;
    end
    
    SelectedFlag(SelectNotYet(idx)) = 1; % select the most informative one
    CountNumRFperImg(imgIdx(SelectNotYet(idx))) = CountNumRFperImg(imgIdx(SelectNotYet(idx))) + 1;
    RFpositive(k) = SelectNotYet(idx);
    
    k = k + 1;
end
fprintf('finished!\n');




