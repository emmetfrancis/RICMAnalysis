function [threshIm] = segmentationRICM(handles)
%segmentationRICM: function to segment RICM images for cellAreaRICM.m
% INPUT: handles structure with (at least) these fields
%   currentImage: image to be segmented 
%   lowVal (: keep pixel values less than this 
%   highVal: keep pixel values higher than this
%   stdVal: keep pixels with a std deviation higher than this in the 3x3
%   neighborhood
% OUTPUT:
% threshIm: final logical segmented image (true = portion of cell)
im = handles.currentImage;
stdIm = stdfilt(im, true(3)); %std deviation matrix for image
try lowVal = str2double(get(handles.lowVal, 'String'));
    highVal = str2double(get(handles.highVal, 'String'));
    stdVal = str2double(get(handles.stdVal, 'String'));
catch
    lowVal = handles.lowVal;
    highVal = handles.highVal;
    stdVal = handles.stdVal;
end
threshIm = im < lowVal | im > highVal; % intensity threshold
threshIm = threshIm | stdIm > stdVal; %std devation threshold
end