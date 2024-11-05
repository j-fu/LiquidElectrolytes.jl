module UnitsConstants
using Unitful
import PhysicalConstants.CODATA2018

sibase(x)=Float64(ustrip(upreferred(x)))
const m=sibase(1Unitful.m)
const K=sibase(1Unitful.K)
const dm=sibase(1Unitful.dm)
const cm=sibase(1Unitful.cm)
const mm=sibase(1Unitful.mm)
const μm=sibase(1Unitful.μm)
const nm=sibase(1Unitful.nm)
const s=sibase(1Unitful.s)
const mA=sibase(1Unitful.mA)
const μA=sibase(1Unitful.μA)
const mol=sibase(1Unitful.mol)
const V=sibase(1Unitful.V)
const mV=sibase(1Unitful.mV)

const AvogadroConstant=sibase(CODATA2018.AvogadroConstant)
const ElementaryCharge=sibase(CODATA2018.ElementaryCharge)
const MolarGasConstant=sibase(CODATA2018.MolarGasConstant)
const FaradayConstant=AvogadroConstant*ElementaryCharge

export m,K,dm,cm,mm,μm,nm,s,mA,μA,mol,V,mV
export AvogadroConstant,ElementaryCharge,MolarGasConstant,FaradayConstant

end
