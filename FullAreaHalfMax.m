function FAHM = FullAreaHalfMax(sourcedistribution,sourcemodel,latency)
%UNTITLED Calculates the full area at half max amplitude
%   Detailed explanation goes here
[~,i_latency] = min(abs(sourcedistribution.time-latency));

% Half max level
half_max = max(sourcedistribution.avg.pow(:,i_latency))/2;

% Find triangles that have at least one point with amplitude >= half max
i_vertices = find(sourcedistribution.avg.pow(:,i_latency)>=half_max);
[triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
triangles = sourcemodel.tri(triangles,:);

% Sum area of triangles and divide by 3 (since its a triangle per point). 
FAHM = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
end