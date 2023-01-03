module StiffnessElements

__precompile__(true)

export stiff_TL1, stiff_QL1

include("legzo.jl")
include("FemElements.jl")
import .FemElements

using LinearAlgebra

## STIFFNESS TL1:
function stiff_TL1(pts, D, thick)
    # Gauss Points (1)
    n_gauss = 1

    x = pts[:, 1]
    y = pts[:, 2]

    # Initialise ke (K element)
    ke = 0

    J::Matrix = zeros(2, 2)
    BB::Matrix = zeros(2, 4)
    B::Matrix = zeros(3, 6)

    # Gauss Points Parameters:
    pg = [1/3]
    w = [0.5]

    # Compute 2D integral of function f
    for i in 1:lastindex(pg)
        xi = pg[i]
        eta = pg[i]
        ww = w[i]

        # Calculate derivative of shape fucntion over (xi,eta)
        dNxi, dNeta = FemElements.TL1_dN(xi, eta)

        # Calculate Jacobian Matrix components using dot() fucntion in order to return ascalar value
        J[1, 1] = dot(dNxi, x)
        J[2, 1] = dot(dNeta, x)
        J[1, 2] = dot(dNxi, y)
        J[2, 2] = dot(dNeta, y)

        # Calculate BB Matrix (Derivative of Shape Functions in cartesian space)
        BB = inv(J) * [dNxi'; dNeta']

        # Assemble B Matrix (Structural Analysis)
        B = [BB[1, 1] 0        BB[2, 1]
             0        BB[2, 1] BB[1, 1]
             BB[1, 2] 0        BB[2, 2]
             0        BB[2, 2] BB[1, 2]
             BB[1, 3] 0        BB[2, 3]
             0        BB[2, 3] BB[1, 3]]

        ke = ke .+ B * D * transpose(B) * thick * det(J) * ww

    end
    return ke
end


## STIFFNESS QL1:
function stiff_QL1(pts, D, thick)
    # Gauss Points (2x2)
    n_gauss = 2

    x = pts[:, 1]
    y = pts[:, 2]

    # Initialise ke (K element)
    ke = 0
    
    J::Matrix = zeros(2,2)
    BB::Matrix = zeros(2,4)
    B::Matrix = zeros(3,8)

    # Gauss Points Calculation:
    pg,w = legzo(n_gauss,1,-1)

    # Compute 2D integral of function f
    for i in 1:lastindex(pg) 
        eta = pg[i]
        w_eta = w[i]
    
        for j in 1:lastindex(pg)
            xi = pg[j]
            w_xi = w[j]

            # Calculate derivative of shape fucntion over (xi,eta)
            dNxi, dNeta = FemElements.QL1_dN(xi,eta)
    
            # Calculate Jacobian Matrix components using dot() fucntion in order to return ascalar value
            J[1,1] = dot(dNxi,x)
            J[2,1] = dot(dNeta,x)
            J[1,2] = dot(dNxi,y)
            J[2,2] = dot(dNeta,y)
            # Calculate BB Matrix (Derivative of Shape Functions in cartesian space)
            BB = inv(J)*[dNxi' ; dNeta']
            
            # Assemble B Matrix (Structural Analysis)
            B = [BB[1,1] 0       BB[2,1] ; 
                 0       BB[2,1] BB[1,1] ;
                 BB[1,2] 0       BB[2,2] ;
                 0       BB[2,2] BB[1,2] ;
                 BB[1,3] 0       BB[2,3] ;
                 0       BB[2,3] BB[1,3] ;
                 BB[1,4] 0       BB[2,4] ;
                 0       BB[2,4] BB[1,4]]

            ke = ke .+ B*D*transpose(B)*thick*det(J)*w_xi*w_eta
        end
    end
    return ke
end

end