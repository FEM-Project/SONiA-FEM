module MatLib

__precompile__(true)

export material

# Material Type:
struct elastMat
    E
    nu
    thick
end

# Elastic Material Library:
function material(MAT_NAME::String)
    if MAT_NAME == "Steel"
        mat = elastMat(210000,0.3,1)
    elseif MAT_NAME == "Alluminum"
        mat = elastMat(70000,0.33,1)
    elseif MAT_NAME == "Custom"
        mat = elastMat(1,0.3,1)
    end
    return mat
end

end