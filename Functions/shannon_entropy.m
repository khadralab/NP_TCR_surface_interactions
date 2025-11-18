function entropy = shannon_entropy(X)
    nbins = max(X) - min(X) + 1;

    counts = histcounts(X, nbins);

    probs = counts ./ sum(counts);
    probs = probs(probs~=0);

    entropy = -sum(probs.*log2(probs));
end