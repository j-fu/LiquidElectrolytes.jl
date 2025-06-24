"""
    AbstractCellData

Abstract type for cell data. Users may implement their own subtypes.

"""
abstract type AbstractCellData end

"""
    electrolytes(::AbstractCellData)

Return a vector of electrolytes, one electrolyte for each cell 
region. 

Subtypes of AbstractCellData can implement their own method for this.
"""
function electrolytes(d::AbstractCellData)
    return d.electrolytes
end

"""
   working_electrode(d::AbstractCellData)

Return boundary index of working electrode.
Subtypes of AbstractCellData can implement their own method for this.
"""
working_electrode(d::AbstractCellData) = working_electrode(electrolytes(d)[1])

"""
   bulk_electrode(d::AbstractCellData)

Return boundary index of bulk electrode.
Subtypes of AbstractCellData can implement their own method for this.
"""
bulk_electrode(d::AbstractCellData) = bulk_electrode(electrolytes(d)[1])

"""
    norm_weights(d::AbstractCellData)

Weights of components in norm calculation.       
Subtypes of AbstractCellData can implement their own method for this.
"""
norm_weights(d::AbstractCellData) = norm_weights(electrolytes(d)[1])

"""
    working_electrode_voltage(d::AbstractCellData)

Return voltage applied at working electrode.
Subtypes of AbstractCellData can implement their own method for this.
"""
working_electrode_voltage(d::AbstractCellData) = working_electrode_voltage(electrolytes(d)[1])

"""
    working_electrode_voltage!(d::AbstractCellData, v)

Set voltage applied at working electrode.
Subtypes of AbstractCellData can implement their own method for this.
"""
working_electrode_voltage!(d::AbstractCellData, v) = working_electrode_voltage!(electrolytes(d)[1], v)

"""
    pressure_index(d::AbstractCellData)

Index of pressure in species list.
Subtypes of AbstractCellData can implement their own method for this.
"""
pressure_index(d::AbstractCellData) = pressure_index(electrolytes(d)[1])

"""
    voltage_index(d::AbstractCellData)

Index of voltage in species list.
Subtypes of AbstractCellData can implement their own method for this.
"""
voltage_index(d::AbstractCellData) = voltage_index(electrolytes(d)[1])

"""
   check_celldata(d::AbstractCellData)

Check if subtypes of AbstractCellData implement all necessary methods.
"""
function check_celldata(d::AbstractCellData)
    return all(
        x -> isa(x, Number),
        [
            working_electrode(d),
            bulk_electrode(d),
            norm_weights(d)[1],
            working_electrode_voltage(d),
            pressure_index(d),
            voltage_index(d),
        ]
    )
end

function update_derived!(d::AbstractCellData)
    check_celldata(d) || error("Incomplete API methods for $(typeof(d))")
    map(update_derived!, electrolytes(d))
    return nothing
end


