function CMoverlap = measures_overlap_prctile(CM1,CM2,percentile)
%MEASURES_OVERLAP_PRCTILE Summary of this function goes here
%   Detailed explanation goes here
totWeights = size(CM1,1)*size(CM1,2);
[rankCM1,sorted_idx1] = sort(abs(CM1(:)),'ascend');
[rankCM2,sorted_idx2] = sort(abs(CM2(:)),'ascend');

rank_idx1 = zeros(1,totWeights);
rank_idx2 = zeros(1,totWeights);
for idx = 1:totWeights
    rank_idx1(idx) = find(sorted_idx1==idx);
    rank_idx2(idx) = find(sorted_idx2==idx);
end

avg_rank = 0.5*(rank_idx1+rank_idx2);

CMoverlap = reshape(avg_rank,size(CM1,1),size(CM1,2));

thresh = prctile(avg_rank(:),percentile);
CMoverlap = (CMoverlap>thresh);

end

