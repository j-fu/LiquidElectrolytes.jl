"""
    AbstractCellData

Abstract type for cell data.
"""
abstract type AbstractCellData end

"""
    electrolytes(::AbstractCellData)

Return a vector of electrolytes, one electrolyte for each cell 
region. Subtypes of AbstractCellData should implement his.
"""
function electrolytes(d::AbstractCellData)
    return d.electrolytes
end

function update_derived!(d::AbstractCellData)
    map(update_derived!, electrolytes(d))
    return nothing
end


working_electrode(d::AbstractCellData) = working_electrode(electrolytes(d)[1])
bulk_electrode(d::AbstractCellData) = bulk_electrode(electrolytes(d)[1])
norm_weights(d::AbstractCellData) = norm_weights(electrolytes(d)[1])
working_electrode_voltage(d::AbstractCellData) = working_electrode_voltage(electrolytes(d)[1])
working_electrode_voltage!(d::AbstractCellData, v) = working_electrode_voltage!(electrolytes(d)[1], v)
pressure_index(d::AbstractCellData) = pressure_index(electrolytes(d)[1])
voltage_index(d::AbstractCellData) = voltage_index(electrolytes(d)[1])
