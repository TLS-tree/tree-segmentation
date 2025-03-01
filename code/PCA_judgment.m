function [rangeX, convexArea] = PCA_judgment(currentTreePoints)

   
    [~, score, ~] = pca(currentTreePoints(:,1:2)); 
   
    rangeX = max(score(:,1)) - min(score(:,1));

    k = convhull(currentTreePoints(:,1), currentTreePoints(:,2));
    convexArea = polyarea(currentTreePoints(k,1), currentTreePoints(k,2));

end