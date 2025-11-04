function [offsets, trl_btn] = findButtonOffsets(trl_in,button_trig,offs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Extract go trials
%tmp = trl(trl(:,4)==13|trl(:,4)==5,:);

%button_trig = raw.trial{1}(contains(raw.label,'STI102'),:);
trig = button_trig>0.5;

trig_smpls = find(trig(2:end)&~trig(1:end-1)) +1;

trl_btn = [];
i_btn = 0;
for i_smpl = 1:length(trig_smpls)
    if any(button_trig(trig_smpls(i_smpl)+(1:10))==button_trig(trig_smpls(i_smpl)))
        i_btn = i_btn+1;
        trl_btn(i_btn,1) = trig_smpls(i_smpl);
    end
end

for i_btn = 1:length(trl_btn)
    [~,i_min] = min(abs((trl_in(:,1)-(trl_in(:,3)))-trl_btn(i_btn)));
    offsets(i_btn,1) = trl_btn(i_btn) - trl_in(i_min,1) + trl_in(i_min,3);
end
offsets = offsets(offsets>0); % remove unmatched or too late presses
offsets = offsets + offs;
trl_btn = trl_in;
trl_btn(:,1:2) = trl_btn(:,1:2) + offsets;
end