% Simply computes euclidean distance
function d = euclidianDistance(x,y)
    d = sqrt(sum((x-y).^2));
end