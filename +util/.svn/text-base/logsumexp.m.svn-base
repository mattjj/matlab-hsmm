function s = logsumexp(a)
if size(a,2) < 2
    s = a;
else
    y = max(a,[],2);
    a = bsxfun(@minus,a,y);
    s = y + log(sum(exp(a),2));
    s(~isfinite(y)) = y(~isfinite(y));
end
end