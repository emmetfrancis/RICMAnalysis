function imOut = polyCrop(imIn, polyPos)
    boundary = boundaryFromVertices(polyPos(:,1), polyPos(:,2));
    x = boundary(:,1);
    y = boundary(:,2);
    % first, determine crop rectangle
    cropRect = [min(x)-1, min(y)-1, max(x)-min(x)+2, max(y)-min(y)+2];
    % now, set all other pixels not in ROI to 'NaN' (appear as black, not
    % considered in intensity or variance analysis
    imLogic = false(size(imIn,1),size(imIn,2));
    for i = 1:length(x)
        imLogic(y(i),x(i)) = true;
    end
    imLogic = imfill(imLogic,'holes');
    imIn(~imLogic) = NaN;
    imOut = imcrop(imIn,cropRect);
end