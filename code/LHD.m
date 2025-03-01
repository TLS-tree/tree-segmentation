function output_LHD = LHD(pointCloudMatrix)

    uniqueLabels = unique(pointCloudMatrix(:, 4));
   
    output_LHD = zeros(length(uniqueLabels), 5);
    
    for i = 1:length(uniqueLabels)
        label = uniqueLabels(i);
       
        currentPoints = pointCloudMatrix(pointCloudMatrix(:, 4) == label, :);
     
        meanX = mean(currentPoints(:, 1));
        meanY = mean(currentPoints(:, 2));
        %meanZ = mean(currentPoints(:, 3));
       
        maxHeightZ = max(currentPoints(:, 3));
        
        distances = sqrt((currentPoints(:, 1) - meanX).^2 + (currentPoints(:, 2) - meanY).^2);
        
        DBH = 2*mean(distances);
       
        output_LHD(i, :) = [label, meanX, meanY, maxHeightZ, DBH ];
    
    end
end