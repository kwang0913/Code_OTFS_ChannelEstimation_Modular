function [Hhat, state] = estimate_orgbl_mimo_2024_v2(Rpilots, PsiPilots, initial, params)
%ESTIMATE_ORGBL_MIMO_2024_V2 Run ORGBL on v2 main-channel pilot dictionaries.
    [Hhat, state] = estimate_orgbl_mimo(Rpilots, PsiPilots, initial, params);
end
