module MatLib

include("InputManager.jl")
import .InputManager
export material

# Material Type:
struct elastMat
    E
    nu
    thick
end

# Elastic Material Library:
function material(MAT_NAME::String)
    E,nu,thick = InputManager.readMat(MAT_NAME)
    mat = elastMat(E,nu,thick)
    return mat
end

end