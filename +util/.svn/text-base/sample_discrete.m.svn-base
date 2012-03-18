function sample = sample_discrete(vec)
    % vec is a column vector to be sampled like a discrete distribution
    % does not need to be normalized
    % vec is interpreted as p(1) = vec(1), p(2) = vec(2), etc.
    cumvals = cumsum(vec,1);
    sample = sum(rand() * cumvals(end) > cumvals) + 1;
end