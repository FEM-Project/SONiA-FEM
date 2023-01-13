module ProblemType

__precompile__(true)

export Plane_Stress, Plane_Strain

function Plane_Stress(young, poiss, thick)
    cc = young / (1 - poiss^2)
    D = cc * [1 poiss 0; poiss 1 0; 0 0 (1-poiss)/2]
    return D
end

function Plane_Strain(young, poiss)
    thick = 1;
    cc = young/((1+poiss)*(1-2*poiss));
    D = cc*[1-poiss poiss 0 ; poiss 1-poiss 0 ; 0 0 (1-2*poiss)/2 ];
    return D
end

end