# Calculate Linear Gauss Points:
function legzo(n::Integer,a,b)
    x::Matrix = zeros(1, n)
    w::Matrix = zeros(1, n)
    m = (n+1)/2
    h = b-a

    for ii::Integer = 1:m
        z = cos(pi*(ii-0.25)/(n+0.5))
        z1 = z+1
        while abs(z-z1)>eps()
            p1 = 1
            p2 = 0
            for jj::Integer = 1:n
                p3 = p2
                p2 = p1
                p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj # The Legendre polynomial.
            end
            global pp = n*(z*p1-p2)/(z^2-1) # The L.P. derivative.
            z1 = z
            z = z1-p1/pp
        end
        x[ii] = z; # Build up the abscissas.
        x[n+1-ii] = -z
        w[ii] = h/((1-z^2)*(pp^2)); # Build up the weights.
        w[n+1-ii] = w[ii]
    end

    if a != -1 || b != 1
        x = (x.+1).*(h/2) .+ a;
    end

    return x,w
end