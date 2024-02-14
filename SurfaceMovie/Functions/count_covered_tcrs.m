function [nt, d] = count_covered_tcrs(tcr_pos, np_pos, rNP)
    d = dist(np_pos', tcr_pos);
    nt = sum(d<rNP,2)';
    d = d<rNP;
end