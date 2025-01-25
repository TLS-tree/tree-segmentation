%% Input trunk slice point cloud after normal vector constraints
lasReader = lasFileReader('F:/tree_segmentation/data/Plot_1.las');
[ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","ScanAngle");
voxelTrunk = ptCloud.Location; 

%% Automatic estimation of Eps and MinPts based on the distance gap of neighborhood point clouds
K = 80;
Mdl = KDTreeSearcher(voxelTrunk);
[minIdxsC,D] = knnsearch(Mdl,voxelTrunk,'K', K);
D1 = D(:,2:K);
Differences = D1(:, 2:end) - D1(:, 1:end-1); 
[maxGaps, indices] = max(Differences, [], 2); 
MeanmaxGaps = mean(maxGaps);
Meanindices = mean(indices+1); 
MinPts = ceil(Meanindices);
eps = MeanmaxGaps*(MinPts-1); 
idx = dbscan(voxelTrunk, eps, MinPts);

[unique_vals, ~, ic] = unique(idx);
counts = accumarray(ic, 1); 
vals_to_replace = unique_vals(counts < 430); 
is_less_than = ismember(idx, vals_to_replace);
idx(is_less_than) = -1; 

ptcloud1 = [voxelTrunk idx];
ptcloud1 = ptcloud1(ptcloud1(:, 4) ~= -1, :);

%% Segmentation result assessment
pointCloudMatrix = double(ptcloud1);

qualifiedPointClouds = [];
unqualifiedPointClouds = [];
treeLabels = unique(pointCloudMatrix(:,4));

for i = 1:length(treeLabels)
    currentTreePoints = pointCloudMatrix(pointCloudMatrix(:,4) == treeLabels(i), 1:3);
    [rangeX, convexArea] = PCA_judgment(currentTreePoints);
    
    % Determine whether the point cloud is valid or invalid based on the conditions
    if rangeX < 1.2 && convexArea > (pi*(1.2/2)^2)/3
        qualifiedPointClouds = [qualifiedPointClouds; currentTreePoints, treeLabels(i)*ones(size(currentTreePoints,1),1)];
    else
        unqualifiedPointClouds = [unqualifiedPointClouds; currentTreePoints, treeLabels(i)*ones(size(currentTreePoints,1),1)];
        
    end
end

%% Perform Hough transform circle detection on incorrect results, and optimize by using the circle center as the new tree position
output_LHD1 = LHD(qualifiedPointClouds);
output_LHD2 = [];

treeLabels2 = unique(unqualifiedPointClouds(:,4));
for i = 1:length(treeLabels2)
    currentTreePoints2 = unqualifiedPointClouds(unqualifiedPointClouds(:,4) == treeLabels2(i), 1:3);
    [centersPointCloud, radii_m] = Hough_Transform(currentTreePoints2, 0.01, 4, 30, 0.88);
    if ~isempty(centersPointCloud)
        centers = [treeLabels2(i)*ones(size(centersPointCloud, 1), 1), centersPointCloud, radii_m];
        %centers = centers(1, :); 
        output_LHD2 = [output_LHD2; centers];
    end    
end

output_LHD = [output_LHD1; output_LHD2];  

