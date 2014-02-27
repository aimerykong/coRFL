function distance = PEDcolor(RFA, RFB, alpha, beta, gamma)
% demo -- show the PED of a pair of receptive fields
%
% ICML2014 paper
%
%
%   Shu Kong (Aimery)
%   aimerykong@gmail.com
%   www.aimerykong.me
%   Dec. 2013

%% prelude
pyramidDesign = 3;

stepA = floor(size(RFA)/pyramidDesign);
stepB= floor(size(RFB)/pyramidDesign);

%% calculate PED of a pair RFs over color feature
distance = zeros( pyramidDesign^2, 1 );
idx = 1;
for i = 1:pyramidDesign
    for j = 1:pyramidDesign
        A = mean(RFA( 1+(i-1)*stepA(1):i*stepA(1), 1+(j-1)*stepA(2):j*stepA(2), 2:3), 1);
        A = mean(A,2);
        B = mean(RFB( 1+(i-1)*stepB(1):i*stepB(1), 1+(j-1)*stepB(2):j*stepB(2), 2:3), 1);
        B = mean(B,2);
        distance(idx) = sum((A(:)-B(:)).^2);
    end
end

distance = sum(sqrt(distance));

%{
%% prelude
pyramidDesign = 3;

%% calculate PED of a pair RFs over color feature
distance = zeros( sum(pyramidDesign), 1 );

spanA = floor(min(size(RFA(:,:,1))) / pyramidDesign(1))-2;
spanB = floor(min(size(RFB(:,:,1))) / pyramidDesign(1))-2;
centerA = floor(size(RFA(:,:,1)) / 2);
centerB = floor(size(RFB(:,:,1)) / 2); 

maskA = zeros(size(RFA,1), size(RFA,2), 2);
tmp = 1:size(RFA,1);
maskA(:,:,1) = repmat(tmp',1,size(RFA,2));
tmp = 1:size(RFA,2);
maskA(:,:,2) = repmat(tmp,size(RFA,1),1);

maskA(:,:,1) = maskA(:,:,1) - centerA(1);
maskA(:,:,2) = maskA(:,:,2) - centerA(2);
maskA = maskA.^2;
maskA = sum(maskA,3);
maskA = sqrt(maskA);

maskB = zeros(size(RFB,1), size(RFB,2), 2);
tmp = 1:size(RFB,1);
maskB(:,:,1) = repmat(tmp',1,size(RFB,2));
tmp = 1:size(RFB,2);
maskB(:,:,2) = repmat(tmp,size(RFB,1),1);

maskB(:,:,1) = maskB(:,:,1) - centerB(1);
maskB(:,:,2) = maskB(:,:,2) - centerB(2);
maskB = maskB.^2;
maskB = sum(maskB,3);
maskB = sqrt(maskB);

for i = 1:pyramidDesign
    innerA = (i-1)*spanA;
    innerB = (i-1)*spanB;
    outerA = i*spanA;
    outerB = i*spanB;
    idxA = find( maskA >= innerA & maskA < outerA);
    idxB = find( maskB >= innerB & maskB < outerB);
    
    tmp = RFA(:,:,1);
    AA(1) = mean(tmp(idxA));
    tmp = RFA(:,:,2);
    AA(2) = mean(tmp(idxA));
    tmp = RFA(:,:,3);
    AA(3) = mean(tmp(idxA));
    
    tmp = RFB(:,:,1);
    BB(1) = mean(tmp(idxB));
    tmp = RFB(:,:,2);
    BB(2) = mean(tmp(idxB));
    tmp = RFB(:,:,3);
    BB(3) = mean(tmp(idxB));
    
    distance(i) = alpha*abs(AA(1)-BB(1)) + ...
        beta*abs(AA(2)-BB(2)) + ...
        gamma*abs(AA(3)-BB(3));
end
    %}
%{
idx = 1;
for pyrIdx = 1:length(pyramidDesign)
    RFsize = size(RFA);
    RFsize = floor( RFsize(1:2) ./ pyramidDesign(pyrIdx) );
    hstepA = RFsize(1);
    wstepA = RFsize(2);
    
    RFsize = size(RFB);
    RFsize = floor( RFsize(1:2) ./ pyramidDesign(pyrIdx) );
    hstepB = RFsize(1);
    wstepB = RFsize(2);
    
    for h = 1:pyramidDesign(pyrIdx)
        for w = 1:pyramidDesign(pyrIdx)
            A = RFA( 1+(h-1)*hstepA:h*hstepA, 1+(w-1)*wstepA:w*wstepA, :);
            B = RFB( 1+(h-1)*hstepB:h*hstepB, 1+(w-1)*wstepB:w*wstepB, :);
            
            A = reshape(A, [size(A,1)*size(A,2), size(A,3)] );
            A = mean(A,1);
            B = reshape(B, [size(B,1)*size(B,2), size(B,3)] );
            B = mean(B,1);
            
            distance(idx) = alpha*norm(A(1)-B(1),'fro') + ...
                beta*norm(A(2)-B(2),'fro') + ...
                gamma*norm(A(3)-B(3),'fro');
            % distance(idx) = min(Graph(:));
            idx = idx + 1;
        end
    end
end
%}
%distance = sum(distance);
