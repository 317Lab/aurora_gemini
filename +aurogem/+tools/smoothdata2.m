%% Copied from later versions of MATLAB

function [B, winsz] = smoothdata2(A, varargin)
% SMOOTHDATA2   Smooth noisy data in two dimensions
%   B = SMOOTHDATA2(A) for a matrix A returns a smoothed version of A using
%   a moving average with a fixed window size.
%
%   B = SMOOTHDATA2(A,METHOD) smooth the entries of A using the specified
%   moving window method METHOD. METHOD must be:
%
%       "movmean"     - (default) smooths by averaging over each window of A.
%       "movmedian"   - smooths by computing the median over each window of A.
%       "gaussian"    - smooths A using a Gaussian filter.
%       "lowess"      - smooths by computing a linear regression over each
%                       window of A.
%       "loess"       - smooths by computing a quadratic regression over
%                       each window of A.
%       "sgolay"      - smooths A using a Savitzky-Golay filter.
%
%   B = SMOOTHDATA2(A,METHOD,WINSIZE) specifies the moving window size used
%   for METHOD. WINSIZE can be a scalar, a two-element cell array of
%   scalars, or a two-element cell array of two-element vectors. By default,
%   WINSIZE is determined automatically from the entries of A.
%
%   B = SMOOTHDATA2(...,NANFLAG) specifies how NaN values are treated:
%
%       "omitmissing" / "omitnan"         -
%                        (default) NaN elements in the input data are ignored
%                        in each window. If all input elements in any window
%                        are NaN, the result for that window is NaN.
%       "includemissing" / "includenan"   -
%                        NaN values in the input data are included when
%                        computing within each window, resulting in NaN.
%
%   B = SMOOTHDATA2(...,SmoothingFactor=FACTOR) specifies a smoothing
%   factor that adjusts the level of smoothing by tuning the default
%   window size. FACTOR must be between 0 (producing smaller moving window
%   lengths and less smoothing) and 1 (producing larger moving window
%   lengths and more smoothing). By default, FACTOR = 0.25. The smoothing
%   factor cannot be specified if WINSIZE input is given.
%
%   B = SMOOTHDATA2(...,SamplePoints={X Y}) also specifies the sample
%   points X and Y used by the smoothing method. X and Y must be numeric,
%   duration, or datetime vectors. They must be sorted and contain unique
%   points. By default, SMOOTHDATA2 uses data sampled uniformly at points
%   X = [1 2 3 ... ] and Y = [1 2 3 ... ]. When 'SamplePoints' are
%   specified, the moving window length is defined relative to the sample
%   points. If X or Y is a duration or datetime vector, then the
%   corresponding moving window size must be a duration.
%
%   B = SMOOTHDATA2(...,'sgolay',...,Degree=D) specifies the degree for
%   the Savitzky-Golay filter. D must be a nonnegative integer.
%
%   [B,WINSIZE] = SMOOTHDATA2(...) also returns the moving window length.
%
%   Example: Smooth a noisy surface
%       A = peaks + randn(size(A))*0.25;
%       B = smoothdata2(A);
%       nexttile
%       surf(A)
%       nexttile
%       surf(B)
%
%   Example: Smooth 2D data using LOESS smoothing
%       A = peaks + randn(size(A))*0.25;
%       B = smoothdata2(A,"loess",7);
%       nexttile
%       surf(A)
%       nexttile
%       surf(B)
%
%   Example: Smooth nonuniform 2D data using Gaussian smoothing and an
%   asymmetric window
%      X = (1:20).^1.2;
%      Y = cumsum(2*sin(pi*(1:20)./21));
%      A = 2*sin(X'/10).*cos(Y/10) + randn(20)*0.25;
%      B = smoothdata2(A,"gaussian",{[4 5],[3 2]},SamplePoints={X,Y});
%      nexttile
%      surf(A)
%      nexttile
%      surf(B)
%
%   See also SMOOTHDATA, FILLMISSING2

%   Copyright 2023 The MathWorks, Inc.

if ~isnumeric(A) && ~islogical(A)
    error(message("MATLAB:smoothdata2:badArray"));
end
if isinteger(A) && ~isreal(A)
    error(message("MATLAB:smoothdata2:complexIntegers"));
end

sparseInput = issparse(A);
if sparseInput
    A = full(A);
end

A = convertToFloat(A);
is2D = true;
[method,winsz,missingFlag,samplePoints,degree] = matlab.internal.math.parseSmoothdataInputs(A, is2D, varargin{:});

% Quick return for empty inputs
if isempty(A)
    B = A;
    if sparseInput
        B = sparse(B);
    end
    return
end

if isempty(samplePoints)
    defaultSamplePoints = true;
    spX = 1:size(A,1);
    spY = 1:size(A,2);
else
    defaultSamplePoints = false;
    spX = samplePoints{1};
    spY = samplePoints{2};
end

dispatchToHelper = startsWith(method,'mov') || matches(method,'gaussian');

spMustBeDouble = false;
spX = matlab.internal.math.convertSamplePoints(spX,spMustBeDouble);
spY = matlab.internal.math.convertSamplePoints(spY,spMustBeDouble);

% Split window sizes
if isscalar(winsz)
    [winX,incUpperBoundX] = matlab.internal.math.convertAndSplitWindow(winsz,spX);
    [winY,incUpperBoundY] = matlab.internal.math.convertAndSplitWindow(winsz,spY);
    if matches(method,"gaussian")
        % Preserve unrounded window size for Gaussian kernel
        if isnumeric(winsz)
            % Use winsz to avoid rounding
            winszX = double(winsz);
            winszY = double(winsz);
        else
            winszX = double(winX);
            winszY = double(winY);
        end
    end
else
    [winX,incUpperBoundX] = matlab.internal.math.convertAndSplitWindow(winsz{1},spX);
    [winY,incUpperBoundY] = matlab.internal.math.convertAndSplitWindow(winsz{2},spY);
    if matches(method,"gaussian")
        % Preserve unrounded window size for calculating Gaussian kernel
        if isnumeric(winsz{1})
            % Use winsz{1} to avoid rounding if winX (cast to sample points type)
            % has lower precision than winsz{1}
            winszX = double(winsz{1});
        else
            winszX = double(winX);
        end

        if isnumeric(winsz{2})
            % Use winsz{2} to avoid rounding if winY (cast to sample points type)
            % has lower precision than winsz{2}
            winszY = double(winsz{2});
        else
            winszY = double(winY);
        end
    end
end

% Send movmean, movmedian, and gaussian methods to helpers
if dispatchToHelper
    if matches(method,"movmean")
        B = matlab.internal.math.movmean2(A,winX,winY,incUpperBoundX,incUpperBoundY,spX,spY,matches(missingFlag,"omitnan"));
    elseif matches(method,"movmedian")
        B = matlab.internal.math.movmedian2(A,winX,winY,incUpperBoundX,incUpperBoundY,spX,spY,matches(missingFlag,"omitnan"));
    else % matches(method,"gaussian")
        if isinteger(spX)
            spX = spX - spX(1); % Avoid artifically large sample points due to conversion
        end
        if isinteger(spY)
            spY = spY - spY(1); % Avoid artifically large sample points due to conversion
        end
        twoSigmaSquaredX = 2*(sum(winszX)/5)^2;
        twoSigmaSquaredY = 2*(sum(winszY)/5)^2;
        B = matlab.internal.math.gaussianSmoothing2(A,winX,winY,incUpperBoundX, ...
            incUpperBoundY,spX,spY,matches(missingFlag,"omitnan"),twoSigmaSquaredX,twoSigmaSquaredY);
    end
else
    if defaultSamplePoints
        spXIsUniform = size(A,1) > 1;
        spYIsUniform = size(A,2) > 1;
    else
        spXIsUniform = isStrictlyUniform(spX);
        spYIsUniform = isStrictlyUniform(spY);
    end
    switch method
        case "loess"
            degree = 2;
            weightedRegression = true;
        case "lowess"
            degree = 1;
            weightedRegression = true;
        case "sgolay"
            % degree passed from parseInputs
            weightedRegression = false;
    end
    B = applyToWindow(A,degree,weightedRegression,missingFlag,winX,winY,incUpperBoundX,incUpperBoundY,spX,spY,spXIsUniform,spYIsUniform);
end
if sparseInput
    B = sparse(B);
end
end

%--------------------------------------------------------------------------

function B = applyToWindow(A,degree,weightedRegression,missingFlag,winX,winY,incUpperBoundX,incUpperBoundY, ...
    spX,spY,spXIsUniform,spYIsUniform)

if spXIsUniform && spYIsUniform && ...
    (matches(missingFlag,"includenan") || ~anynan(A))
    xSpacing = spX(2)-spX(1);
    kx = floor(winX/xSpacing);
    ySpacing = spY(2)-spY(1);
    ky = floor(winY/ySpacing);

    % Account for exclusive upper bounds
    if ~incUpperBoundX && mod(winX(2),xSpacing) == 0
        kx(2) = max(0,kx(2)-1);
    end
    if ~incUpperBoundY && mod(winY(2),ySpacing) == 0
        ky(2) = max(0,ky(2)-1);
    end

    if kx(1)+kx(2)+1 < size(A,1) && ky(1)+ky(2)+1 < size(A,2)
        % Use the convolution code path if the window is smaller than the
        % data
        B = localRegressionUniformSP(A,kx,ky,xSpacing,ySpacing,degree,weightedRegression);
        return;
    end
    % If the window is larger than the data, flow through to the general
    % code path
end

% Precompute windows
[xWindowStarts, xWindowEnds] = generateWindowBounds(winX,spX,incUpperBoundX);
[yWindowStarts, yWindowEnds] = generateWindowBounds(winY,spY,incUpperBoundY);

[numRows,numCols,numPages] = size(A);
% Preallocate B
B = A;
for jj = 1:numCols
    for ii = 1:numRows
        for kk = 1:numPages
            xStart = xWindowStarts(ii);
            xEnd   = xWindowEnds(ii);
            yStart = yWindowStarts(jj);
            yEnd   = yWindowEnds(jj);
            B(ii,jj,kk) = localRegression(A(xStart:xEnd,yStart:yEnd,kk),[xStart,xEnd],[yStart,yEnd], ...
                ii,jj,spX,spY,degree,weightedRegression,missingFlag);
        end
    end
end
end

%--------------------------------------------------------------------------

function [windowStarts, windowEnds] = generateWindowBounds(window,sp,incUpperBound)
windowStarts = zeros(1,numel(sp));
windowEnds = zeros(1,numel(sp));
numSp = numel(sp);
% Find last window that includes the first element
ind = 1;
while (ind + 1 <= numSp) && (sp(1) + window(1) >= sp(ind + 1))
    ind = ind + 1;
end
windowStarts(1:ind) = 1;
% Find the last element in the first window
lastElement = ind;
currentWindowSPEnd = sp(ind) + window(2);
if incUpperBound
    while (lastElement + 1 <= numSp) && (currentWindowSPEnd >= sp(lastElement + 1))
        % Increment lastElement until it is outside the window
        lastElement = lastElement + 1;
    end
else
    while (lastElement + 1 <= numSp) && (currentWindowSPEnd > sp(lastElement + 1))
        % Increment lastElement until it is outside the window
        lastElement = lastElement + 1;
    end
end
windowEnds(1:ind) = lastElement;
ind = ind + 1;

firstElement = 1;
while ind <= numSp && lastElement < numSp
    while sp(firstElement) + window(1) < sp(ind)
        firstElement = firstElement + 1;
    end
    windowStarts(ind) = firstElement;

    lastElement = max(lastElement,ind);
    currentWindowSPEnd = sp(ind) + window(2);
    if ~incUpperBound
        if isinteger(sp)
            % This is safe for unsigned types because window(2) cannot be
            % zero in this branch, so currentWindowSPEnd is >= 1
            currentWindowSPEnd = currentWindowSPEnd - 1;
        else
            currentWindowSPEnd = currentWindowSPEnd - eps(currentWindowSPEnd);
        end
    end
    while (lastElement + 1 <= numSp) && (currentWindowSPEnd >= sp(lastElement + 1))
        % Increment lastElement until it is outside the window
        lastElement = lastElement + 1;
    end
    windowEnds(ind) = lastElement;
    ind = ind + 1;
end
windowStarts(ind:end) = firstElement;
windowEnds(min(ind,numSp):end) = numSp; % Guarantees that the last sample point is always in the last window
end

%--------------------------------------------------------------------------

function b = localRegression(data,xWindowBounds,yWindowBounds,ii,jj,spX,spY,degree,useWeights, ...
    missingFlag)
% Convert numeric sample points to doubles to avoid integer negative
% overflow
if isnumeric(spX)
    spX = double(spX);
end
if isnumeric(spY)
    spY = double(spY);
end

omitNaN = missingFlag == "omitnan";
xDist = getRelativeDistance(xWindowBounds,ii,spX);
yDist = getRelativeDistance(yWindowBounds,jj,spY);

if useWeights
    weights = getTricubicWeights(xDist,yDist);
else
    weights = [];
end

if omitNaN
    nanLocations = isnan(data);
    if all(nanLocations)
        b = nan;
        return
    end
    dataPointsToInclude = ~nanLocations(:);
else
    dataPointsToInclude = true(numel(data),1);
end

if useWeights
    dataPointsToIncludeTmp = dataPointsToInclude & weights~=0;
    weights = weights(dataPointsToIncludeTmp);
    if numel(weights) == 0
        useWeights = false;
    else
        dataPointsToInclude = dataPointsToIncludeTmp;
    end
end
data = data(dataPointsToInclude);
data = data(:);

V = buildVandermondeMatrix(xDist,yDist,degree,useWeights,weights,dataPointsToInclude);

if useWeights
    weightedData = double(data).*weights(:);
else
    weightedData = double(data);
end
bTmp = matlab.internal.math.leastSquaresFit(V,weightedData);
b = bTmp(1);
end

%--------------------------------------------------------------------------

function b = localRegressionUniformSP(data,kx,ky,xSpacing,ySpacing,degree,useWeights)
% Reshape data to (at most) 3 dimensions so that the output shape from
% convn is at most 3D
dataShape = size(data);
data = data(:,:,:);

hasNaNs = anynan(data);
% Remove NaNs before convolving
% conv2 does not ignore NaNs when the corresponding weight is 0
if hasNaNs && useWeights
    nanLocations = isnan(data);
    data(nanLocations) = zeros('like',data);
    replaceWithNaN = false(size(data));
end

winszX = kx(1)+kx(2)+1;
winszY = ky(1)+ky(2)+1;

b = data;
b = localRegressionCenterBlock(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b);
if hasNaNs && useWeights
    replaceWithNaN = localRegressionCenterBlock(nanLocations,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,replaceWithNaN);
end

symmetricEqualWindow = kx(1) == kx(2) && isequal(kx,ky) && xSpacing == ySpacing;
if symmetricEqualWindow
    k = 2*kx(1)+1;
    b = localRegressionSymmetricEqualCorners(data,k,degree,useWeights,b);
    b = localRegressionSymmetricEqualSides(data,k,degree,useWeights,b);
    if hasNaNs && useWeights
        replaceWithNaN = localRegressionSymmetricEqualCorners(nanLocations,k,degree,useWeights,replaceWithNaN);
        replaceWithNaN = localRegressionSymmetricEqualSides(nanLocations,k,degree,useWeights,replaceWithNaN);
        b(replaceWithNaN) = NaN;
    end
else
    b = localRegressionXEdges(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b);
    b = localRegressionYEdges(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b);
    b = localRegressionCorners(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b);
    if hasNaNs && useWeights
        replaceWithNaN = localRegressionXEdges(nanLocations,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,replaceWithNaN);
        replaceWithNaN = localRegressionYEdges(nanLocations,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,replaceWithNaN);
        replaceWithNaN = localRegressionCorners(nanLocations,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,replaceWithNaN);
        b(replaceWithNaN) = NaN;
    end
end
b = reshape(b,dataShape);
end

%--------------------------------------------------------------------------

function b = localRegressionXEdges(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b)
[xSize,ySize] = size(data,[1 2]);
if kx > xSize
    % There are no locations where the filter is invariant in the
    % x-direction
    return
end
ii = kx(1) + 1; % This is a dummy x-coordinate in the middle of the block of data
sp = {(1:winszX)*xSpacing,(1:winszY)*ySpacing};
ky = min(ky,ySize);
% Left edge
for jj = 1:ky(1)
    h = buildLocalRegressionFilter(winszX,winszY,ii,jj,sp,degree,useWeights);
    b(kx(1)+1:xSize-kx(2),jj,:) = convn(data(:,1:winszY,:),h,'valid');
end

% Right edge
for jj = ky(1)+2:winszY
    h = buildLocalRegressionFilter(winszX,winszY,ii,jj,sp,degree,useWeights);
    b(kx(1)+1:xSize-kx(2),ySize-winszY+jj,:) = convn(data(:,ySize-winszY+1:end,:),h,'valid');
end
end

%--------------------------------------------------------------------------

function b = localRegressionYEdges(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b)
[xSize,ySize] = size(data,[1 2]);
if ky > ySize
    % There are no locations where the filter is invariant in the
    % y-direction
    return
end
jj = ky(1) + 1; % This is a dummy x-coordinate in the middle of the block of data
sp = {(1:winszX)*xSpacing,(1:winszY)*ySpacing};
kx = min(kx,xSize);
% Top edge
for ii = 1:kx(1)
    h = buildLocalRegressionFilter(winszX,winszY,ii,jj,sp,degree,useWeights);
    b(ii,ky(1)+1:ySize-ky(2),:) = convn(data(1:winszX,:,:),h,'valid');
end

% Bottom edge
for ii = kx(1)+2:winszX
    h = buildLocalRegressionFilter(winszX,winszY,ii,jj,sp,degree,useWeights);
    b(xSize-winszX+ii,ky(1)+1:ySize-ky(2),:) = convn(data(xSize-winszX+1:end,:,:),h,'valid');
end
end

%--------------------------------------------------------------------------

function b = localRegressionCorners(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b)
[xSize,ySize] = size(data,[1 2]);
spX = (1:winszX)*xSpacing;
spY = (1:winszY)*ySpacing;
windowBoundsX = [1 winszX];
windowBoundsY = [1 winszY];
kx = min(kx,xSize);
ky = min(ky,ySize);

for kk = 1:size(data,3)
    topLeftData = data(1:winszX,1:winszY,kk);
    topRightData = data(1:winszX,ySize-winszY+1:end,kk);
    for ii = 1:kx(1)
        % Top left corner
        for jj = 1:ky(1)
            b(ii,jj,kk) = localRegression(topLeftData,windowBoundsX,windowBoundsY,ii,jj,spX,spY,degree, ...
                useWeights,'includeNaN');
        end

        % Top right corner
        for jj = ky(1)+2:winszY
            b(ii,ySize-winszY+jj,kk) = localRegression(topRightData,windowBoundsX,windowBoundsY,ii,jj, ...
                spX,spY,degree,useWeights,'includeNaN');
        end
    end

    bottomLeftData = data(xSize-winszX+1:end,1:winszY,kk);
    bottomRightData = data(xSize-winszX+1:end,ySize-winszY+1:end,kk);
    for ii = kx(1)+2:winszX
        % Bottom left corner
        for jj = 1:ky(1)
            b(xSize-winszX+ii,jj,kk) = localRegression(bottomLeftData,windowBoundsX,windowBoundsY,ii,jj, ...
                spX,spY,degree,useWeights,'includeNaN');
        end

        % Bottom right corner
        for jj = ky(1)+2:winszY
            b(xSize-winszX+ii,ySize-winszY+jj,kk) = localRegression(bottomRightData,windowBoundsX, ...
                windowBoundsY,ii,jj,spX,spY,degree,useWeights,'includeNaN');
        end
    end
end
end

%--------------------------------------------------------------------------

function b = localRegressionSymmetricEqualCorners(data,k,degree,useWeights,b)
sizeOfBlock = (k-1)/2; % k must be odd if the window is symmetric
[xSize,ySize,numPages] = size(data);
sp = {1:min(k,xSize),1:min(k,ySize)};
for ii = 1:sizeOfBlock
    for jj = ii:sizeOfBlock
        h = buildLocalRegressionFilter(k,k,ii,jj,sp,degree,useWeights);
        h = flip(flip(h,1),2);
        h = h(:)';
        for kk = 1:numPages
            b(ii,jj,kk)                         = h*              reshape(data(1:k,   1:k,                 kk),[],1); % Top left corner
            b(xSize - ii + 1,jj,kk)             = h*reshape(data(xSize:-1:(xSize-k+1),1:k,                 kk),[],1); % Bottom left corner
            b(ii,ySize - jj + 1,kk)             = h*reshape(data(1:k,                 ySize:-1:(ySize-k+1),kk),[],1); % Top right corner
            b(xSize - ii + 1,ySize - jj + 1,kk) = h*reshape(data(xSize:-1:(xSize-k+1),ySize:-1:(ySize-k+1),kk),[],1); % Bottom right corner
            if ii ~= jj
                hT = reshape(reshape(h,k,k)',1,[]);
                b(jj,ii,kk)                         = hT*reshape(data(1:k,                 1:k,                 kk),[],1); % Top left corner
                b(xSize - jj + 1,ii,kk)             = hT*reshape(data(xSize:-1:(xSize-k+1),1:k,                 kk),[],1); % Bottom left corner
                b(jj,ySize - ii + 1,kk)             = hT*reshape(data(1:k,                 ySize:-1:(ySize-k+1),kk),[],1); % Top right corner
                b(xSize - jj + 1,ySize - ii + 1,kk) = hT*reshape(data(xSize:-1:(xSize-k+1),ySize:-1:(ySize-k+1),kk),[],1); % Bottom right corner
            end
        end
    end
end
end

%--------------------------------------------------------------------------

function b = localRegressionSymmetricEqualSides(data,k,degree,useWeights,b)
sizeOfBlock = (k-1)/2; % k must be odd if the window is symmetric
[xSize,ySize] = size(data,[1 2]);
jj = sizeOfBlock + 1; % This is a dummy y-coordinate in the middle of the block of data
sp = {1:min(k,xSize),1:min(k,ySize)};
for ii = 1:sizeOfBlock
    h = buildLocalRegressionFilter(k,k,ii,jj,sp,degree,useWeights);
    b(ii,(sizeOfBlock + 1):(ySize-sizeOfBlock),:)             = convn(data(1:k,:,:),                h,      'valid'); % Top
    b((sizeOfBlock + 1):(xSize-sizeOfBlock),ii,:)             = convn(data(:,1:k,:),                h',     'valid'); % Left
    b(xSize - ii + 1,(sizeOfBlock + 1):(ySize-sizeOfBlock),:) = convn(data((xSize - k + 1):end,:,:),flip(h),'valid'); % Bottom
    b((sizeOfBlock + 1):(xSize-sizeOfBlock),ySize - ii + 1,:) = convn(data(:,(ySize - k + 1):end,:),flip(h)','valid'); % Right
end
end

%--------------------------------------------------------------------------

function b = localRegressionCenterBlock(data,kx,ky,xSpacing,ySpacing,degree,useWeights,winszX,winszY,b)
[xSize,ySize] = size(data,[1 2]);
ii = kx(1) + 1; % This is a dummy x-coordinate in the middle of the block of data
jj = ky(1) + 1; % This is a dummy y-coordinate in the middle of the block of data
sp = {(1:winszX)*xSpacing,(1:winszY)*ySpacing};
h = buildLocalRegressionFilter(winszX,winszY,ii,jj,sp,degree,useWeights);
b(kx(1) + 1:xSize - kx(2), ky(1) + 1:ySize - ky(2),:) = convn(data,h,'valid');
end

%--------------------------------------------------------------------------

function d = getRelativeDistance(windowBounds,queryIndex,sp)
% Convert numeric sample points to doubles to avoid integer negative
% overflow
sp = double(sp);
d = sp(windowBounds(1):windowBounds(2)) - sp(queryIndex);
end

%--------------------------------------------------------------------------

function weights = getTricubicWeights(xDist,yDist)
if numel(xDist) == 1 && numel(yDist) == 1
    % Quick return when the is only one point in the window
    weights = 1;
    return
end
totalDist = sqrt((xDist(:).^2+yDist.^2));
totalDist = totalDist./max(totalDist,[],'all');
weights = sqrt(((1 - totalDist(:).^3).^3));

if all(weights == 0)
    weights = ones(numel(totalDist),1);
end
end

%--------------------------------------------------------------------------

function V = buildVandermondeMatrix(xDist,yDist,degree,useWeights,weights,dataPointsToInclude)
% Check that there are enough elements in each direction to perform the
% requested regression (this is a necessary but not sufficient conditon)
reshapedPointsToInclude = reshape(dataPointsToInclude,[numel(xDist),numel(yDist)]);
xDegree = min(degree, max(sum(any(reshapedPointsToInclude,2),1)) - 1);
yDegree = min(degree, max(sum(any(reshapedPointsToInclude,1),2)) - 1);

xDegree = max(xDegree,0);
yDegree = max(yDegree,0);

Vx = power(reshape(xDist,[],1),0:xDegree);
Vy = power(reshape(yDist,[],1),0:yDegree);
V = kron(Vy,Vx);
V = V(dataPointsToInclude,:);
if useWeights
    V = V.*weights;
end

% Remove columns associated with polynomial degrees that exceed the
% specified degree
columnPolyDegree = reshape((0:yDegree)'+(0:xDegree),1,[]);
columnsToDrop = columnPolyDegree > degree;
V(:,columnsToDrop) = [];
end

%--------------------------------------------------------------------------

function h = buildLocalRegressionFilter(sizeX,sizeY,ii,jj,sp,degree,useWeights)
xDist = getRelativeDistance([1 sizeX],ii,sp{1});
yDist = getRelativeDistance([1 sizeY],jj,sp{2});
numDataPoints = numel(xDist)*numel(yDist);

if useWeights
    weights = getTricubicWeights(xDist,yDist);
    dataPointsToInclude = weights~=0;
    weights(~dataPointsToInclude) = [];
else
    weights = [];
    dataPointsToInclude = true(1,numDataPoints);
end

V = buildVandermondeMatrix(xDist,yDist,degree,useWeights,weights,dataPointsToInclude);

[Q,~] = qr(V,'econ');
queryIndInWin = (jj-1)*sizeX + ii;
queryIndInQ = sum(dataPointsToInclude(1:queryIndInWin)); % Compensate for omitted data points
hInternal = Q*Q(queryIndInQ,:)';
hInternal = hInternal(:);

if useWeights
    hInternal = hInternal.*weights;
end
h = zeros(1,numDataPoints);
h(dataPointsToInclude) = hInternal;
h = reshape(h,numel(xDist),numel(yDist));
h = flip(flip(h,1),2);
end

%--------------------------------------------------------------------------

function A = convertToFloat(A)
if isobject(A)
    if isequal(underlyingType(A),'single')
        A = single(A);
    else
        A = double(A);
    end
elseif isinteger(A)
    A = double(A);
end
end

%--------------------------------------------------------------------------

function tf = isStrictlyUniform(sp)
if numel(sp) < 2
    tf = false;
    return
end
diffSp = diff(sp); % Sample points must be in ascending order, so this is safe for unsigned types
tf = all(diffSp == diffSp(1));
end