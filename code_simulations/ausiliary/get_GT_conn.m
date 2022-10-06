function [gtConn,gtDelay] = get_GT_conn(N1,N2,D,post11,post12,post21,post22,delays11,delays12,delays21,delays22)
%% Ground truth connectivity matrix
gtConn = zeros(N1+N2,N1+N2);
gtDelay = nan(N1+N2,N1+N2);

% Within M1 conn
gtConn11 = zeros(N1,N1);
gtDelay11 = nan(N1,N1);
for sIdx = 1:N1 %sender
    for rIdx = 1:N1 %receiver
        tmpConn = sum(post11(sIdx,:) == rIdx);
        if tmpConn
            gtConn11(sIdx, rIdx) = tmpConn;
        end
    end
    for dIdx = 1:D
        tmpDel = delays11(sIdx,dIdx);
        if ~isempty(tmpDel{1})
            for nIdx = 1:numel(tmpDel{1})
                gtDelay11(sIdx, post11(sIdx,tmpDel{1}(nIdx))) = dIdx;
            end
        end
    end
end

% Within DLS conn
gtConn22 = zeros(N2,N2);
gtDelay22 = nan(N2,N2);
for sIdx = 1:N2 %sender
    for rIdx = 1:N2 %receiver
        tmpConn = sum(post22(sIdx,:) == rIdx);
        if tmpConn
            gtConn22(sIdx, rIdx) = tmpConn;
        end
    end
    for dIdx = 1:D
        tmpDel = delays22(sIdx,dIdx);
        if ~isempty(tmpDel{1})
            for nIdx = 1:numel(tmpDel{1})
                gtDelay22(sIdx, post22(sIdx,tmpDel{1}(nIdx))) = dIdx;
            end
        end
    end
end

% M1 to DLS conn
gtConn12 = zeros(N1,N2);
gtDelay12 = nan(N1,N2);
for sIdx = 1:N1 %sender
    for rIdx = 1:N2 %receiver
        tmpConn = sum(post11(sIdx,:) == rIdx);
        if tmpConn
            gtConn12(sIdx, rIdx) = tmpConn;
        end
    end
    for dIdx = 1:D
        tmpDel = delays12(sIdx,dIdx);
        if ~isempty(tmpDel{1})
            for nIdx = 1:numel(tmpDel{1})
                gtDelay12(sIdx, post12(sIdx,tmpDel{1}(nIdx))) = dIdx;
            end
        end
    end
end

% M1 to DLS conn
gtConn21 = zeros(N2,N1);
gtDelay21 = nan(N2,N1);
for sIdx = 1:N2 %sender
    for rIdx = 1:N1 %receiver
        tmpConn = sum(post11(sIdx,:) == rIdx);
        if tmpConn
            gtConn21(sIdx, rIdx) = tmpConn;
        end
    end
    for dIdx = 1:D
        tmpDel = delays21(sIdx,dIdx);
        if ~isempty(tmpDel{1})
            for nIdx = 1:numel(tmpDel{1})
                gtDelay21(sIdx, post21(sIdx,tmpDel{1}(nIdx))) = dIdx;
            end
        end
    end
end

gtConnFrom1 = [gtConn11;gtConn21];
gtConnFrom2 = [gtConn12;gtConn22];
gtDelFrom1 = [gtDelay11;gtDelay21];
gtDelFrom2 = [gtDelay12;gtDelay22];

gtConn = [gtConnFrom1,gtConnFrom2];
gtDelay = [gtDelFrom1,gtDelFrom2];
end

