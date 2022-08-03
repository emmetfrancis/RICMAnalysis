function varargout = tracesFromBinary(binIm, minSize, blur)
%tracesFromBinary: find continuous outline in a binary image
%   INPUT:
%   binIm: input image as a binary matrix
%   OUTPUT:
%   varargout{1}: n-by-2 matrix, where column 1 is row values, column 2 are
%   column values, and n is the number of points to outline the image
%   varargout{2}: information from regionprops after filling holes and
%   smoothing
numRow = size(binIm, 1);
numCol = size(binIm, 2);
if any(binIm(:))
    % first, blur
    if blur > 0
        binIm = imfill(binIm, 'holes');
        range = -round(blur):round(blur);
        [xMesh, yMesh] = meshgrid(range,range);
%         kernel = blur - sqrt(xMesh.^2 + yMesh.^2);
        kernel = exp(-(xMesh.^2+yMesh.^2)/(2*(blur/3)^2));
        kernel(kernel<0) = 0;
        kernel = kernel ./ sum(kernel(:));
        binIm = im2double(binIm);
        binIm = conv2(binIm, kernel, 'same');
        binIm = imbinarize(binIm);
    end
    % fill in holes and choose largest area
    binIm = imfill(binIm, 'holes'); 
%     stats = regionprops(binIm, 'Area', 'PixelIdxList');
    connComps = bwconncomp(binIm);
    PixelLists = connComps.PixelIdxList;
    areas = zeros(1,length(PixelLists));
    for k = 1:length(PixelLists)
        areas(k) = length(PixelLists{k}); % num pixels in region
    end
%     areas = [stats(:).Area];
    n = length(areas);
    [~, idx] = max(areas);
    excl = horzcat(1:(idx-1), idx+1:n); % indices for shapes to be filled with black
    for i = 1:length(excl)
%         binIm(stats(excl(i)).PixelIdxList) = false;
        binIm(PixelLists{excl(i)}) = false;
    end
    if max(areas) < minSize
        binIm = false(size(binIm));
    end
    if nargout == 2
        varargout{2} = binIm;
    end
    % now, locate edge pixel
    foundTrue = false;
    for row = 1:numRow
        for col = 1:numCol
            if binIm(row,col)
                initialRow = row;
                initialCol = col;
                foundTrue = true;
                break
            end
        end
        if foundTrue
            break
        end
    end
    if foundTrue
        initialPt = [initialRow, initialCol];
        % trace varargout{1}
%         varargout{1} = bwtraceboundary(binIm, initialPt, 'E');
        allBound = bwboundaries(binIm);
        if length(allBound) > 1
            error('Other nonzero regions not exluded correctly')
        end
        varargout{1} = allBound{1};
    else % likely because shape not larger than minSize
        varargout{1} = [0, 0];
    end
else
    % no true pixels -> no outline
    varargout{1} = [0, 0];
end
if nargout == 2
    varargout{2} = binIm;
end
end