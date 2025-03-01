%% Input trunk position and normalized point cloud
root = load('F:/.../data/output_LHD.csv'); 
lasReader = lasFileReader('F:/.../data/Plot_subset.las'); 
[ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","ScanAngle");
p1 = ptCloud.Location;

%% l0 Cut Pursuit supervoxel segmentation
PR_REG_STRENGTH1=0.1;     
PR_MIN_NN1=3;              

PR_EDGE_STRENGTH1=1;       
PR_DECIMATE_RES1=0.05;     
option=1;
speed=0;
verbose=0;
PR_DECIMATE_RES2=0.1;      

k = 10;                    
w_max = 0.2;               % Similarity weight threshold

segs=[p1(:,1:end),zeros(size(p1,1),1)];      
subPos=p1; 
subPos=subPos-mean(subPos,1); 
[~,ia,ic]=unique(floor((subPos-min(subPos,[],1))/PR_DECIMATE_RES1)+1,'rows');
decPos=subPos(ia,:);  
nNodes=size(decPos,1);
if nNodes>0
    inComponent=cutpursuit(decPos,PR_MIN_NN1,PR_REG_STRENGTH1,PR_EDGE_STRENGTH1,option,speed,verbose);
    segs(:,end)=double(inComponent(ic))+1;
end
labels = segs(:, 4);
uniqueLabels = unique(labels);
groupSizes = accumarray(labels, 1); 

sums_x = accumarray(labels, segs(:, 1), [], @sum);
sums_y = accumarray(labels, segs(:, 2), [], @sum);
sums_z = accumarray(labels, segs(:, 3), [], @sum);

clusterCenters1 = zeros(size(uniqueLabels,1),3);
for i = 1:size(uniqueLabels,1)
    label = uniqueLabels(i);
    clusterCenters1(i, :) = [sums_x(label), sums_y(label), sums_z(label)] ./ groupSizes(label);
end
voxelCenters = clusterCenters1;
ptCloudid = segs;

% Add tree trunk position after point cloud voxelization
Ptcloud = [root(:, 2:4); voxelCenters]; 
Ptcloud = double(Ptcloud);
Roots = (1 : size(root, 1)); 
ptCloudid1 = ptCloudid;
ptCloudid1(:, 4) = ptCloudid1(:, 4) + size(root, 1); 
root1 = [root, (1:size(root, 1))'];
ptCloudid1 = [root1(:, [2,3,4,6]); ptCloudid1]; 

%% Construction of edges and weight calculation in the supervoxel graph
segs = ptCloudid1;
segs=[segs,zeros(size(segs,1),1)]; 
segsMin=min(segs(:,1:3),[],1); 
segs(:,1:3)=segs(:,1:3)-segsMin; 

clusterIdx=segs(:,end-1); 
[clusterU,clusterUIdx,clusterV]=unique(clusterIdx); 
[clusterVSorted,clusterVIdx]=sort(clusterV);
clusterVGroup = mat2cell(clusterVIdx, diff([0;find([(diff(clusterVSorted) ~= 0);1])]),1);

currentClusterDecs=zeros(size(segs,1),4)-1; 
currentClusterIs={}; 

% Extract the coordinates of each initial cluster center point
clusterCentroids=zeros(numel(clusterU),3);
iter=0;
for i=1:numel(clusterU)  
    currentClusterPos=segs(clusterVGroup{i},1:3);
    [~,ia,ic]=unique(floor((currentClusterPos-min(currentClusterPos,[],1))/PR_DECIMATE_RES2)+1,'rows');

    startI=iter+1;
    iter=iter+numel(ia);
    end_i=iter;
    currentClusterDecs(startI:end_i,:)=[currentClusterPos(ia,:),repmat(i,numel(ia),1)]; 
    currentClusterIs{i}=startI:end_i;
    clusterCentroids(i,:)=mean(currentClusterPos,1);
end
currentClusterDecs(currentClusterDecs(:,1)==-1,:)=[]; 
currentClusterDecsMean=mean(currentClusterDecs(:,1:3),1); 
currentClusterDecs(:,1:3)=currentClusterDecs(:,1:3)-currentClusterDecsMean; 

Mdl = KDTreeSearcher(clusterCentroids);
[minIdxsC,~] = knnsearch(Mdl,clusterCentroids,'K',k+1);

nnDists=zeros(size(minIdxsC))+10; 

if isempty(gcp('nocreate'))
    parpool;
end

DD = zeros(size(minIdxsC, 1), 3); 
parfor i=1:size(minIdxsC,1)
    currentClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,1)},1:3); 
    currentClusterDec_center = mean(currentClusterDec, 1);
    min_dis1 = sqrt(sum((currentClusterDecs(1:size(root, 1),1:2) - currentClusterDec_center(:,1:2)).^2, 2)); 
    min_dis1 = min(min_dis1);
    for j=2:size(minIdxsC,2)
        nnClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,j)},1:3); 
        nnClusterDec_center = mean(nnClusterDec,1);
        min_dis2 = sqrt(sum((currentClusterDecs(1:size(root, 1),1:2) - nnClusterDec_center(:,1:2)).^2, 2)); 
        min_dis2 = min(min_dis2);
        min_dis = max(min_dis1, min_dis2);

        D_xy = norm(currentClusterDec_center(1:2)-nnClusterDec_center(1:2));
        D_z = abs(currentClusterDec_center(3)-nnClusterDec_center(3));
        D = [D_xy, D_z, min_dis];
        DD(i, :) = D;  
    end
end
max_value = max(DD);

L = size(minIdxsC,2);
parfor i=1:size(minIdxsC,1)
    currentClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,1)},1:3); 
    for j=2:L
        nnClusterDec=currentClusterDecs(currentClusterIs{minIdxsC(i,j)},1:3); 
        w = exp(-(DD(i,1)/max_value(1,1))^2) * exp(-(DD(i,2)/max_value(1,2))^2) * exp(-(DD(i,3)/max_value(1,3))^2);
        if w > w_max
            [~,D]=jitknnsearch(nnClusterDec,currentClusterDec); 
            nnDists(i,j)=min(D);
        else
            nnDists(i,j)=NaN;
        end  
        
    end
end

nnDists(:,1)=0;
nnDists = nnDists(:,2:end);

minIdxsC = minIdxsC(:,2:end); 
n_point = size(Ptcloud,1);
source1 = reshape(repmat(1:n_point, [k 1]), [1 (k * n_point)])';
target1 = reshape(minIdxsC', [1 (k * n_point)])';
Edge = [source1, target1];
weight1 = reshape(nnDists', [], 1);

Edge = sort(Edge,2); 
[Edge, ia, ~] = unique(Edge, 'rows', 'first'); 
weight1 = weight1(ia); 

ia = find(~isnan(weight1)); 
weight1 = weight1(ia);  
Edge = Edge(ia,:);    

Edge2 = [Edge(:,2),Edge(:,1)];
Edge = [Edge;Edge2];
weight2 = weight1;
weight1 = [weight1; weight2];

weight3 = sqrt(sum((Ptcloud(Edge(:,1),:) - Ptcloud(Edge(:,2),:)).^2,2));
weight1 = 0.5*weight1 +0.5*weight3;      % Weighted average of the supervoxel center distance and boundary distance

%% Graph construction
adj = sparse(Edge(:,1),Edge(:,2),weight1,length(Ptcloud),length(Ptcloud));
Gp = graph(adj);

%% shortest path to each root
d = distances(Gp, Roots, 'Method', 'positive');
d(d == inf) = nan;
[~,treeid] = min(d,[],1,'omitnan');
treeid = treeid'; 
Ptcloudid = [Ptcloud, treeid]; 
Ptcloudid = Ptcloudid((size(root, 1)+1):end, :); 
clear adj Gp d DBH treeid

new_column = Ptcloudid(ptCloudid(:,4), 4);  
ptCloudID = [ptCloudid new_column]; % Final segmentation result

XYZ = ptCloudID(:, 1:3);
label = ptCloudID(:, 5);
ptCloudID1 = pointCloud(XYZ); 
ptCloudID1.Intensity = label; % Here, 'intensity' actually represents the point cloud's label information
pcwrite(ptCloudID1, 'F:/.../data/seg_result_subset.ply', 'PLYFormat', 'binary'); 