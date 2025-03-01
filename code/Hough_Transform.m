function [centersPointCloud, radii_m] = Hough_Transform(points, gridSize, rmin, rmax, Sensitivity)

% etermine the image size
xMin = floor(min(points(:,1)));
xMax = ceil(max(points(:,1)));
yMin = floor(min(points(:,2)));
yMax = ceil(max(points(:,2)));
imageWidth = ceil((xMax - xMin) / gridSize);
imageHeight = ceil((yMax - yMin) / gridSize);

binaryImage = zeros(imageHeight, imageWidth);

for i = 1:size(points, 1)
    % Calculate the position of points in the grid
    col = floor((points(i,1) - xMin) / gridSize) + 1;
    row = floor((points(i,2) - yMin) / gridSize) + 1;
    
    % Set the corresponding grid cell value to 1
    binaryImage(row, col) = 1;
end

[centers, radii] = imfindcircles(binaryImage, [rmin rmax], 'ObjectPolarity', 'bright', 'Sensitivity', Sensitivity, 'Method', 'PhaseCode');
radii_m = radii * gridSize;
if ~isempty(centers)
    
    % Convert circle center coordinates to the original point cloud coordinate system
    centersPointCloud = zeros(size(centers, 1), 3);
    for i = 1:size(centers, 1)
    centersPointCloud(i,1) = (centers(i,1) - 1) * gridSize + xMin;
    centersPointCloud(i,2) = (centers(i,2) - 1) * gridSize + yMin;
    centersPointCloud(i,3) = min(points(:,3)); 
    end

else

    centersPointCloud = centers;

end

end

