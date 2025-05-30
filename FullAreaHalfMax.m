function [peak_activation] = FullAreaHalfMax(sourcedistribution,sourcemodel,peak,params,save_path)
%UNTITLED Calculates the full area at half max amplitude
%   Detailed explanation goes here

if length(peak.peak_latency) == 2
    [~,i1] = min(abs(sourcedistribution.time-peak.peak_latency(1)));
    [~,i2] = min(abs(sourcedistribution.time-peak.peak_latency(2)));
    
    dat = sourcedistribution.avg.mom(:,i1:i2).^2;
    [~,i_latency] = max(max(dat,[],1)); % max of max across sources
    i_latency = i1-1+i_latency;
    latency = sourcedistribution.time(i_latency);
else
    [~,i_latency] = min(abs(sourcedistribution.time-peak.peak_latency));
    latency = sourcedistribution.time(i_latency);
end

% Half max level
[peak_mom, i_maxsource] = max(abs(sourcedistribution.avg.mom(:,i_latency)));
half_max = peak_mom/2;
peak_loc = sourcemodel.pos(i_maxsource,:);
peak_pow = sourcedistribution.avg.mom(i_maxsource,i_latency).^2; % power at max latency and source

% Find triangles that have at least one point with amplitude >= half max
i_halfmax_vertices = find(abs(sourcedistribution.avg.mom(:,i_latency))>=half_max);
halfmax_distribution = abs(sourcedistribution.avg.mom(:,i_latency))>=half_max;
[triangles,~] = find(ismember(sourcemodel.tri,i_halfmax_vertices)); 
triangles = sourcemodel.tri(triangles,:);

% Sum area of triangles and divide by 3 (since its a triangle per point).
fahm = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;  

peak_activation = []; 
peak_activation.label = peak.label;
peak_activation.latency = latency;
peak_activation.loc = peak_loc;
peak_activation.pow = peak_pow;
peak_activation.mom = peak_mom;
peak_activation.fahm = fahm;
peak_activation.halfmax_distribution = halfmax_distribution;
end