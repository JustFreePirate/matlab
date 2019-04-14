function res = vectorGetLerp( v, idx )
    n = length(v);
    l = floor(idx);
    r = ceil(idx);
    res = lerp(v(l), v(r), idx - l);
end

