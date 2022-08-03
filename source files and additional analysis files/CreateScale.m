function [scaleVal] = CreateScale()
%CreateScale allows the user to select an image of the scale and resize a
%rectangle to define vertical and horizontal dimensions to set the micron
%to pixel ratio
%   No input required, prompts the user to select an image file to use for
%   the function
%   Output: scaleVal is the micron to pixel ratio
[file,path] = uigetfile({'*.bmp';'*.tiff'},'Choose scale image');
if file~=0
    figure()
    scalePic = imread(fullfile(path,file));
    imshow(scalePic)
    rectHandle = imrect;
    waitforbuttonpress
    val=get(gcf,'CurrentKey');
    while ~strcmp(val,'return')
        waitforbuttonpress
        val=get(gcf,'CurrentKey');
    end
    api = iptgetapi(rectHandle);
    position = api.getPosition;
    dimensions = position(3:4);
    micronValCell = inputdlg({'Horizontal Microns','Vertical Microns'},'Micron Values',1,{'N/A','N/A'});
    micronVals = zeros(2,1);
    for i=1:2
        num = str2double(micronValCell{i});
        if ~isnan(num)
            micronVals(i) = num;
        else
            micronVals(i) = 0;
        end
    end
    indices = find(micronVals);
    micronVals = micronVals(indices);
    dimensions = dimensions(indices);
    scaleVal = mean(micronVals ./ dimensions);
    close(gcf)
end
end