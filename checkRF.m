function mat = checkRF(threshold, mat, imgA, imgB, RFlist_A, RFlist_B)
% This function checks whether the fed receptive field (RF) is a valid
% one, and return true if valid or false if invalid.
%
%
%
%
% Shu Kong (Aimery)
% aimerykong@gmail.com
% Feb. 2014

upperbound = 20;
%% intra-contrast of color 
for i = 1:size(RFlist_A,2) % from imageA
    tmpRF = imgA(RFlist_A(1,i):RFlist_A(3,i), RFlist_A(2,i):RFlist_A(4,i));
    %check the current RF
    if sum(tmpRF(:))/numel(tmpRF) < threshold
        mat(i,:) = upperbound;
        mat(:,i) = upperbound;
    end
end

for i = 1:size(RFlist_B,2) % from imageA
    tmpRF = imgB(RFlist_B(1,i):RFlist_B(3,i), RFlist_B(2,i):RFlist_B(4,i));
    %check the current RF
    if sum(tmpRF(:))/numel(tmpRF) < threshold
        mat(i,:) = upperbound;
        mat(:,i) = upperbound;
    end
end
