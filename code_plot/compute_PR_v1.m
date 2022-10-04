function PR_measures = compute_PR_v1(gtConn,peakTE,peakHOTE,peakXCov,peakXCorr,minutes_range,params)

for rIdx = 1:params.nRep
    rIdx
    for mIdx = 1:numel(minutes_range)
        mIdx
        
        idxT = (gtConn{rIdx,mIdx} == 1); %idxs of true connections
        idxF = (gtConn{rIdx,mIdx} == 0); %idxs of false connections
        % CM = connectivity measures
        for idx = 1:99
            if params.doOverlap
                tmpConnOverlap = measures_overlap_prctile(peakXCov{rIdx,mIdx},peakHOTE{rIdx,mIdx},idx);
                PR_measures.HOTE_XCov.tp(rIdx,mIdx,idx) = mean(tmpConnOverlap(idxT)); % TPR
                PR_measures.HOTE_XCov.fp(rIdx,mIdx,idx) = mean(tmpConnOverlap(idxF)); % FPR
                PR_measures.HOTE_XCov.prec(rIdx,mIdx,idx) = sum(sum(tmpConnOverlap(idxT)))/sum(sum(tmpConnOverlap)); % precision = TP/(TP+FP) = TP/(inferred P)
                PR_measures.HOTE_XCov.rec(rIdx,mIdx,idx) = sum(sum(tmpConnOverlap(idxT)))/sum(sum(gtConn{rIdx,mIdx})); % recall = TP/(TP+FN) = TP/(Actual P)
            end
            
            % Compute PR for TE
            TEthr = prctile(peakTE{rIdx,mIdx}(:),idx);
            tmpTEmap = (peakTE{rIdx,mIdx} > TEthr);
            PR_measures.TE.tp(rIdx,mIdx,idx) = mean(tmpTEmap(idxT)); % TPR
            PR_measures.TE.fp(rIdx,mIdx,idx) = mean(tmpTEmap(idxF)); % FPR
            PR_measures.TE.prec(rIdx,mIdx,idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(tmpTEmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            PR_measures.TE.rec(rIdx,mIdx,idx) = sum(sum(tmpTEmap(idxT)))/sum(sum(gtConn{rIdx,mIdx})); % recall = TP/(TP+FN) = TP/(Actual P)
            
            % Compute PR for HOTE
            HOTEthr = prctile(peakHOTE{rIdx,mIdx}(:),idx);
            tmpHOTEmap = (peakHOTE{rIdx,mIdx} > HOTEthr);
            PR_measures.HOTE.tp(rIdx,mIdx,idx) = mean(tmpHOTEmap(idxT)); % TPR
            PR_measures.HOTE.fp(rIdx,mIdx,idx) = mean(tmpHOTEmap(idxF)); % FPR
            PR_measures.HOTE.prec(rIdx,mIdx,idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(tmpHOTEmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            PR_measures.HOTE.rec(rIdx,mIdx,idx) = sum(sum(tmpHOTEmap(idxT)))/sum(sum(gtConn{rIdx,mIdx})); % recall = TP/(TP+FN) = TP/(Actual P)
            
            % Compute PR for XCov
            XCovthr = prctile(abs(peakXCov{rIdx,mIdx}(:)),idx);
            tmpXCovmap = (abs(peakXCov{rIdx,mIdx}) > XCovthr);
            PR_measures.XCov.tp(rIdx,mIdx,idx) = mean(tmpXCovmap(idxT)); % TPR
            PR_measures.XCov.fp(rIdx,mIdx,idx) = mean(tmpXCovmap(idxF)); % FPR
            PR_measures.XCov.prec(rIdx,mIdx,idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(tmpXCovmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            PR_measures.XCov.rec(rIdx,mIdx,idx) = sum(sum(tmpXCovmap(idxT)))/sum(sum(gtConn{rIdx,mIdx})); % recall = TP/(TP+FN) = TP/(Actual P)
            
            % Compute PR for XCov
            XCorrthr = prctile(abs(peakXCorr{rIdx,mIdx}(:)),idx);
            tmpXCorrmap = (abs(peakXCorr{rIdx,mIdx}) > XCorrthr);
            PR_measures.XCorr.tp(rIdx,mIdx,idx) = mean(tmpXCorrmap(idxT)); % TPR
            PR_measures.XCorr.fp(rIdx,mIdx,idx) = mean(tmpXCorrmap(idxF)); % FPR
            PR_measures.XCorr.prec(rIdx,mIdx,idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(tmpXCorrmap)); % precision = TP/(TP+FP) = TP/(inferred P)
            PR_measures.XCorr.rec(rIdx,mIdx,idx) = sum(sum(tmpXCorrmap(idxT)))/sum(sum(gtConn{rIdx,mIdx})); % recall = TP/(TP+FN) = TP/(Actual P)
        end
    end
end