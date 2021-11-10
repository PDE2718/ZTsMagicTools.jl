module ZTsMagicTools

export axshpy!, PLv, QLv

function axshpy!(a, x, y; biasx = 0)
    n = length(x)
    b = mod(biasx,n)
    @inbounds @simd for i=1:n-b
        y[i+b] += a*x[i]
    end
    @inbounds @simd for i=1:b
        y[i] += a*x[n-b+i]
    end
end

@inline function PLv(x::T, l::Int) where T
    P = zeros(T,l+1)
    P[1]=1. ; P[2]=x
    for n=3:l+1
        @inbounds P[n]=((2n-3)*x*P[n-1]-(n-2)*P[n-2])/(n-1)
    end
    return P
end

@inline function QLv(x::T, l::Int) where T
    Q = zeros(T,l+1)
    Q[1]=log((1+x)/(1-x))/2; Q[2]=x*Q[1]-1
    for n=3:l+1
        @inbounds Q[n]=((2n-3)*x*Q[n-1]-(n-2)*Q[n-2])/(n-1)
    end
    return Q
end

end # module
