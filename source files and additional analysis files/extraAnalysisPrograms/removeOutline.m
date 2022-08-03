[file, path] = uigetfile('*.mat', 'Choose Cell Prints File');
loadStruct = load(fullfile(path,file));
try cellPrints = loadStruct.cellPrints;
catch cellPrints = loadStruct.CellPrints;
end
cellPrints = ~cellPrints;
for i = 1:size(cellPrints,3)
    binIm = cellPrints(:,:,i);
    outline = bwperim(binIm);
    binIm(outline) = false;
    cellPrints(:,:,i) = binIm;
end
cellPrints = ~cellPrints;
newFile = strcat('RemovedOL', file);
save(fullfile(path,newFile), 'cellPrints')