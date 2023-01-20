module MatLib

include("PreProcessing.jl")
import .PreProcessing
export material

# Material Type:
struct elastMat
    E
    nu
    thick
end

# Elastic Material Library:
function material(MAT_NAME::String)
    E,nu,thick = PreProcessing.readMat(MAT_NAME)
    mat = elastMat(E,nu,thick)
    return mat
end

end