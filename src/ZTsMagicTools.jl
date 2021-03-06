module ZTsMagicTools

export axshpy!, PLv, QLv, fx2y!, fxy2z!

function axshpy!(a, x, y; biasx = 0)
    n = length(x)
    b = mod(biasx,n)
    @inbounds begin
        for i=1:n-b
            y[i+b] += a*x[i]
        end
        for i=1:b
            y[i] += a*x[n-b+i]
        end
    end
end

function PLv(l::Int, x::T) where T
    P = zeros(T,l+1)
    P[1]=1. ; P[2]=x
    for n=3:l+1
        @inbounds P[n]=((2n-3)*x*P[n-1]-(n-2)*P[n-2])/(n-1)
    end
    return P
end

function QLv(l::Int, x::T) where T
    Q = zeros(T,l+1)
    Q[1]=log((1+x)/(1-x))/2; Q[2]=x*Q[1]-1
    for n=3:l+1
        @inbounds Q[n]=((2n-3)*x*Q[n-1]-(n-2)*Q[n-2])/(n-1)
    end
    return Q
end

function fx2y!(f, x, y)                     #  y = f(x)
    if size(x)==size(y)
        @inbounds for i = 1:length(x)
            y[i] = f(x[i])
        end
    else
        throw(ErrorException("Dimension not match"))
    end
end

function fxy2z!(f, x, y, z)                 # z = f(x,y)
    if size(x)==size(y)==size(z)
        @inbounds for i = 1:length(x)
            z[i] = f(x[i],y[i])
        end
    else
        throw(ErrorException("Dimension not match"))
    end
end

end # module
