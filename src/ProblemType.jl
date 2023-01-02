module ProblemType

__precompile__(true)

export Plane_Stress

function Plane_Stress(young, poiss, thick)
    cc = young / (1 - poiss^2)
    D = cc * [1 poiss 0; poiss 1 0; 0 0 (1-poiss)/2]
    return D
end

end