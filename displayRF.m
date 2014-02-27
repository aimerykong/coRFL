function displayRF(SelectedFlag, imgIdx, subPath, RFtemplates)
% visualize the selected RFs

if nargin < 4
    RFtemplates = genRFcandidates();
end
numRFperImg = size(RFtemplates, 2);

IdxImg = imgIdx(SelectedFlag); % image index
IdxImg = IdxImg(:);
IdxRF = SelectedFlag(:) - (IdxImg-1)*numRFperImg;   % RF index in this image
idxSet = [IdxImg IdxRF];

DIRpath = dir(subPath);
disp(idxSet);
for i = size(idxSet, 1):-1:1
    im = imread( fullfile(subPath, DIRpath(idxSet(i)+2).name) );
    im = im2double(im);
    imSize = size(im);
    if max(imSize) > 400
        im = imresize(im, 400/max(imSize)); % resize image with preserved aspect ratio
        imSize = size(im);
    end
    
    RFi = RFtemplates(1:4, idxSet(i, 2));
    RFi([1 3], :) = floor( (RFi([1 3], :)-1)*(imSize(1)/100))+1;
    RFi([2 4], :) = floor( (RFi([2 4], :)-1)*(imSize(2)/100))+1;
    
    RFplus = im(RFi(1):RFi(3), RFi(2):RFi(4), :);
    im( :, :, 1) = 1;%im( :, :, 1)*0.5;
    im( :, :, 2) = im( :, :, 2)*0.8;
    im(:, :, 3) = im( :, :, 3)*0;
    im(RFi(1):RFi(3), RFi(2):RFi(4), :) = RFplus; %im(RFi(1):RFi(3), RFi(2):RFi(4), :) .* mask;
    clear mask;
    
    figure; imshow(im); title(['selected RF-' num2str(i)]);
end
