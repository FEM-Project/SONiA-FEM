module FemElements

__precompile__(true)

export EM1_N, EM1_dN, EM2_N, EM2_dN, TL1_N, TL1_dN, QL1_N, QL1_dN

using LinearAlgebra

## EM1:
# Shape Functions
function EM1_N(xi)
    N = zeros(2)
    N[1] = (1-xi)/2
    N[2] = (1+xi)/2
    return N
end

# Derivative of Shape Functions
function EM1_dN(xi)
    dN = zeros(2)
    dN[1] = -1/2
    dN[2] = 1/2
    return dN
end


## EM2 :
# Shape Functions
function EM2_N(xi)
    N = zeros(3)
    N[1] = xi/2*(xi-1)
    N[2] = (1+xi)*(1-xi)
    N[3] = xi/2*(xi+1)
    return N
end

# Derivative of Shape Functions
function EM2_dN(xi)
    dN = zeros(3)
    dN[1] = xi-1/2
    dN[2] = -2*xi
    dN[3] = xi+1/2
    return dN
end


## TL1:
# Shape Function
function TL1_N(xi, eta)

    # Domain definition:
    p1 = [0 0]
    p2 = [1 0]
    p3 = [0 1]

    n = [   1 p1[1] p1[2] ;
            1 p2[1] p2[2] ;
            1 p3[1] p3[2]  ]

    if rank(n) < size(n)[1]
        println("[ERROR] Non compatible element!")
    else
        N = [1 xi eta]/n
    end

    return N[:]
end

# Derivative of Shape Function
function TL1_dN(xi,eta)
    
    # Domain definition:
    p1 = [0 0]
    p2 = [1 0]
    p3 = [0 1]

    n = [   1 p1[1] p1[2] ;
            1 p2[1] p2[2] ;
            1 p3[1] p3[2]  ]

    if rank(n) < size(n)[1]
        println("[ERROR] Non compatible element!")
    else
        N = inv(n)
    end

    dNxi = transpose(N[2,:])
    dNeta = transpose(N[3,:])

    return dNxi[:], dNeta[:]

end


## QL1:
# Shape Function
function QL1_N(xi,eta)

    # Domain definition:
    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]

    n = [   1 p1[1] p1[2] p1[1]*p1[2] ;
            1 p2[1] p2[2] p2[1]*p2[2] ;
            1 p3[1] p3[2] p3[1]*p3[2] ;
            1 p4[1] p4[2] p4[1]*p4[2]   ]

    if rank(n) < size(n)[1]
        println("[ERROR] Non compatible element!")
    else
        N = [1 xi eta xi*eta]/n
    end

    return N[:]

end

# Derivative of Shape Function
function QL1_dN(xi,eta)
    
    # Domain definition:
    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]

    n = [   1 p1[1] p1[2] p1[1]*p1[2] ;
            1 p2[1] p2[2] p2[1]*p2[2] ;
            1 p3[1] p3[2] p3[1]*p3[2] ;
            1 p4[1] p4[2] p4[1]*p4[2]   ]

    if rank(n) < size(n)[1]
        println("[ERROR] Non compatible element!")
    else
        N = inv(n)
    end

    Nxi = transpose([N[2,:] N[4,:]])
    dNxi = [1 eta]*Nxi

    Neta = transpose([N[3,:] N[4,:]])
    dNeta = [1 xi]*Neta

    return dNxi[:], dNeta[:]

end

end
