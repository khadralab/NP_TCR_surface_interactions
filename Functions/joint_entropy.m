function je = joint_entropy(X,Y)

    nbins_X = max(X) - min(X) + 1;
    nbins_Y = max(Y) - min(Y) + 1;

    counts = histcounts2(X,Y, [nbins_X, nbins_Y]);

    probs = counts ./ sum(counts, "all");
    probs = probs(probs~=0);

    je = -sum(probs.*log2(probs));
end