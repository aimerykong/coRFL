function [RFlist imgCenter] = genRFcandidates( type, LeapStepFactor, imgSize, PartitionFactor )
%   Given an image and a decisive factor, the function generates a set of
%   recpetive field candidates, and records the partition of these
%   regions/fields in a string.
%
% input
%   imgSize                 - the pixel resolution of image (size);
%   PartitionFactor     - overcomplete rectangular bins as grids that
%                                   constitutes a receptive field candidate; 
%
% output
%   RFlist                     - a 4-by-N matrix that index all the N receptive field
%                                   candidates; the 4 integers of each
%                                   column record positions of upper-left
%                                   and lower-right corners;
%                                   
%
%   Shu Kong (Aimery)
%   aimerykong@gmail.com
%   www.aimerykong.me
%   Nov. 2013

%% default parameters to create a receptive field (RF) model 
if nargin < 4
    PartitionFactor = 4;
end
if nargin < 3    
    imgSize = [100 100];
end
if nargin < 2    
    LeapStepFactor = 0.5;
end
if nargin < 1
    type = 'pyramid'; % 'grid'; %
end

imgSize = imgSize(1:2);
imgCenter = round(imgSize/2);

%%
gridSize = floor(imgSize / PartitionFactor);
RFlist = [];

if strcmp(type,'pyramid')    
    binSize = PartitionFactor; % small grid size
    LeapStep = floor( gridSize*LeapStepFactor );
    for ii = 1:binSize
        for jj = 1: binSize
            regionSize = gridSize + ([ii jj]-1).*gridSize;
            
            for vertical = 1:PartitionFactor/LeapStepFactor-ii+1
                anchor(1) = 1 + LeapStep(1)*(vertical-1);
                if anchor(1) > imgSize(1)
                    continue;
                end
                for horizontal = 1:PartitionFactor/LeapStepFactor-jj+1
                    anchor(2) = 1 + LeapStep(2)*(horizontal-1);
                    
                    posBin = zeros(4,1);
                    posBin(1:2) = anchor(:);
                    posBin(3:4) = anchor(:) + regionSize(:)-1;
                    posBin(5:6) = round( (posBin(1:2)+posBin(3:4)) /2 );
                    posBin(7) = norm(abs(posBin(5:6) - imgCenter(:)), 'fro');
                    if posBin(3) <= imgSize(1) && posBin(4) <= imgSize(2)
                        RFlist = [RFlist, posBin];
                    end
                    
                end
            end
        end
    end
elseif strcmp(type,'grid')
    binSize = floor(PartitionFactor/2); % small grid size
    LeapStep = floor( gridSize*LeapStepFactor );
    for ii = 1:binSize
        for jj = 1: binSize
            regionSize = gridSize + ([ii jj]-1).*gridSize;
            
            for vertical = 1:PartitionFactor/LeapStepFactor-ii+1
                anchor(1) = 1 + LeapStep(1)*(vertical-1);
                if anchor(1) > imgSize(1)
                    continue;
                end
                for horizontal = 1:PartitionFactor/LeapStepFactor-jj+1
                    anchor(2) = 1 + LeapStep(2)*(horizontal-1);
                    
                    posBin = zeros(4,1);
                    posBin(1:2) = anchor(:); % upper-left corner
                    posBin(3:4) = anchor(:) + regionSize(:)-1; % bottom-right corner
                    posBin(5:6) = round( (posBin(1:2)+posBin(3:4)) /2 ); % RF center
                    posBin(7) = norm(abs(posBin(5:6) - imgCenter(:)), 'fro'); % distance to image center
                    if posBin(3) <= imgSize(1) && posBin(4) <= imgSize(2)
                        RFlist = [RFlist, posBin];
                    end                    
                end
            end
        end
    end
else
end



