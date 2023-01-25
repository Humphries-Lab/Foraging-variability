function [SubjMinNLL, SubjMinNLLParams, SubjNLLEval] = GridSearchForaging(SubjData, Env, AlphaVector, BetaVector)

% create grid on which to initialise searches
[AlphaPatchGrid,AlphaRhoGrid, BetaGrid] = ndgrid(AlphaVector,AlphaVector,BetaVector);

nAlpha = size(AlphaPatchGrid, 1);
nBeta = size(AlphaPatchGrid, 3);

SubjNLLEval = zeros(nAlpha,nAlpha,nBeta);

for iAlphaPatch = 1:nAlpha
    for iAlphaRho = 1:nAlpha
        for iBeta = 1:nBeta
            params = [AlphaPatchGrid(iAlphaPatch,iAlphaRho,iBeta),AlphaRhoGrid(iAlphaPatch,iAlphaRho,iBeta), BetaGrid(iAlphaPatch,iAlphaRho,iBeta)];  % get the current set of initial search parameters
            
            [SubjNLLEval(iAlphaPatch,iAlphaRho,iBeta)] = NLL_M1_MVT_RWRho(params, Env, SubjData);                       
            %[SubjNLLEval(iAlphaPatch,iAlphaRho,iBeta)] = NLL_M2_MVT_RW(params, Env, SubjData);            
            %[SubjNLLEval(iAlphaPatch,iAlphaRho,iBeta)] = NLL_M5_RLOn(params, Env, SubjData);            
            %[SubjNLLEval(iAlphaPatch,iAlphaRho,iBeta)] = NLL_M21_RLOff(params, Env, SubjData);
            %[SubjNLLEval(iAlphaPatch,iAlphaRho,iBeta)] = NLL_M25_RLOn(params, Env, SubjData);

            SubjMinNLL = min(min(min(SubjNLLEval)));   % minimum negative log likelihood over all starting positions
            [ixAlphaPatch,ixAlphaRho,ixBeta] = ind2sub(size(AlphaPatchGrid),find(SubjNLLEval == SubjMinNLL));    % indices of location of minimum
            ixAlphaPatch = ixAlphaPatch(1); %Note: if more than one starting location converges on same parameters then will have same NLL, so use first
            ixAlphaRho = ixAlphaRho(1);
            ixBeta = ixBeta(1);

            % get the corresponding fitted parameter values.
            SubjMinNLLParams = [AlphaPatchGrid(ixAlphaPatch,ixAlphaRho,ixBeta), AlphaRhoGrid(ixAlphaPatch, ixAlphaRho,ixBeta),BetaGrid(ixAlphaPatch, ixAlphaRho,ixBeta)];
        end
    end
end
