function boundary = boundaryFromVertices(xpos, ypos)
% this function defines pixel values for a boundary given the vertices in
% pixel values for a polygonal region of interest

xpos = round(xpos);
ypos = round(ypos);
x = xpos(1);
y = ypos(1);
for i = 2:length(xpos)+1
    x1 = xpos(i-1);
    y1 = ypos(i-1);
    if i == length(xpos)+1
        x2 = xpos(1);
        y2 = ypos(1);
    else
        x2 = xpos(i);
        y2 = ypos(i);
    end
    dist = sqrt((x2-x1)^2 + (y2-y1)^2);
    for j = 1:round(dist)-1
        xCur = x1 + j*(x2-x1)/dist;
        yCur = y1 + j*(y2-y1)/dist;
        x = vertcat(x, round(xCur));
        y = vertcat(y, round(yCur));
    end
    if i < length(xpos)+1
        x = vertcat(x, x2);
        y = vertcat(y, y2);
    end
end
boundary = horzcat(x,y);

end

