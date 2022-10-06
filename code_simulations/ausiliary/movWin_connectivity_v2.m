function [TE_time,HOTE_time,XCov_time,XCorr_time] = ...
    movWin_connectivity_v2(firings,TE_pairs,TE_delays,...
                            doHO,HOTE_pairs,HOTE_delays,...
                            doXCov,XCov_pairs,XCov_delays,max_delay,...
                            doXCorr,XCorr_pairs,XCorr_delays,...
                            eps_del,newBin,timeJump,del_consistency,norm_method)
% This function estimates TE over time using a moving window of size newBin
% with jumps of length timeJump. The computation is taken only for pairs
% that resulted significant under the full-time TE analysis and for the
% peak delay (CONSIDER adding with a precision eps over delays)
% out: TE_time = #pairs x #downsampled time points

% YOU HAVE TO RESTRUCTURE THE INPUT TO THIS FUNCTION grouping "do"s, "pairs" and "delays" into structures 

nUnits = size(firings,1);
oldTimePoints = size(firings,2);

newTimePoints = floor((oldTimePoints-(newBin-timeJump))/timeJump);

TE_time=zeros(size(TE_pairs,1),newTimePoints); 
HOTE_time=TE_time; 
XCov_time=TE_time; 
XCorr_time=TE_time;

oldTcount=1;
for t = 1:newTimePoints
    t
    disp([num2str(t),' out of ',num2str(newTimePoints)])
    
    tidx = oldTcount:oldTcount+newBin-1;
    for pairIdx = 1:size(TE_pairs,1)
        tmpFirings = firings([TE_pairs(pairIdx,1), TE_pairs(pairIdx,2)],tidx);
        asdf = SparseToASDF(tmpFirings, 1);
        if del_consistency
            delay_range = (TE_delays(pairIdx)-eps_del):(TE_delays(pairIdx)+eps_del);
            delay_range(delay_range<1)=[];
        else
            delay_range = 1:max_delay;
        end
        tmp = ASDFTE(asdf, delay_range, 1, 1); 
        if TE_pairs(pairIdx,1) < TE_pairs(pairIdx,2)
            TE_time(pairIdx,t) = tmp(1,2);
        elseif TE_pairs(pairIdx,1) > TE_pairs(pairIdx,2)
            TE_time(pairIdx,t) = tmp(2,1);
        end   
    end
    % Compute HOTE
    if doHO
        for pairIdx = 1:size(HOTE_pairs,1)
            tmpFirings = firings([HOTE_pairs(pairIdx,1), HOTE_pairs(pairIdx,2)],tidx);
            asdf = SparseToASDF(tmpFirings, 1);
            if del_consistency
                delay_range = (HOTE_delays(pairIdx)-eps_del):(HOTE_delays(pairIdx)+eps_del);
                delay_range(delay_range<1)=[];
            else
                delay_range = 1:max_delay;
            end
            tmp = ASDFTE(asdf, delay_range, 5, 5);
            if HOTE_pairs(pairIdx,1) < HOTE_pairs(pairIdx,2)
                HOTE_time(pairIdx,t) = tmp(1,2);
            elseif HOTE_pairs(pairIdx,1) > HOTE_pairs(pairIdx,2)
                HOTE_time(pairIdx,t) = tmp(2,1);
            end 
        end
    end
    % Compute XCov
    if doXCov
        for pairIdx = 1:size(XCov_pairs,1)
            cell1_st = firings(XCov_pairs(pairIdx,1),tidx);
            cell2_st = firings(XCov_pairs(pairIdx,2),tidx);
            %tmpSpikeTrains = [spikeTrains(i,:);spikeTrains(j,:)]';
            [xCov,lags] = xcov(cell1_st,cell2_st,max_delay,norm_method);
            xCov(ceil(numel(xCov)/2))=[]; % eliminating 0 delays
            lags(ceil(numel(xCov)/2))=[];
            if del_consistency
                selDel=find(lags==(-XCov_delays(pairIdx)));
                %xCov = reshape(xCov,2*max_delay+1,2,2); % elems ordered as [11,21;12,22]
                % Neurons can be either correlated or anticorrelated --> take abs value
                XCov_time(pairIdx,t)=xCov(selDel);
            else
                negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                %xCov = reshape(xCov,2*max_delay+1,2,2); % elems ordered as [11,21;12,22]
                % Neurons can be either correlated or anticorrelated --> take abs value
                [~,tmpxCovDelay]=max(abs(xCov(negIdxs)));
                XCov_time(pairIdx,t)=xCov(tmpxCovDelay);
            end
        end
    end
    if doXCorr
        for pairIdx = 1:size(XCorr_pairs,1)
            cell1_st = firings(XCorr_pairs(pairIdx,1),tidx);
            cell2_st = firings(XCorr_pairs(pairIdx,2),tidx);
            %tmpSpikeTrains = [spikeTrains(i,:);spikeTrains(j,:)]';
            [xCorr,lags] = xcorr(cell1_st,cell2_st,max_delay,norm_method);
            xCorr(ceil(numel(xCorr)/2))=[]; % eliminating 0 delays
            lags(ceil(numel(xCorr)/2))=[];
            if del_consistency
                selDel=find(lags==(-XCorr_delays(pairIdx)));
                %xCov = reshape(xCov,2*max_delay+1,2,2); % elems ordered as [11,21;12,22]
                % Neurons can be either correlated or anticorrelated --> take abs value
                XCorr_time(pairIdx,t)=abs(xCorr(selDel));
            else
                negIdxs = (lags<0);  % computing xcov(a,b) shifts b by lag, so if we find maximum at lag=(-del) it means that the link is a->b with lag = del
                %xCov = reshape(xCov,2*max_delay+1,2,2); % elems ordered as [11,21;12,22]
                % Neurons can be either correlated or anticorrelated --> take abs value
                [~,tmpxCovDelay]=max(abs(xCorr(negIdxs)));
                XCorr_time(pairIdx,t)=xCorr(tmpxCovDelay);
            end
        end
    end
    oldTcount = oldTcount + timeJump;
end

end

