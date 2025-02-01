"""
    Specifications of RRDEs 
"""

mutable struct DiscSpec
    r_1
    r_2
    r_3
    r_4
    w_gap
    h_cyl
    DiscSpec(disktype) = DiscSpec(new(), disktype)
end

function DiscSpec(self, disktype)
    @local_unitfactors μm mm cm
    println("Loading specs for RRDE type $(disktype)")
    if disktype == "E6R1"
        # Pine E6 Series ChangeDisk RRDE Tips
        # https://www.pineinst.com/echem/viewproduct.asp?ID=45651
        # Collection Efficiency: 25%
        # Exposed Tip Length (t): 25.4 mm
        # Disk Insert Thickness (m): 4.0 mm
        self.r_1 = 2.5 * mm  # Disk Outer Diameter (D1): 5.0 mm
        self.r_2 = 3.25 * mm # Ring Inner Diameter (D2): 6.5 mm
        self.r_3 = 3.75 * mm # Ring Outer Diameter (D3): 7.5 mm
        self.r_4 = 5 * mm   # Shroud Outer Diameter (s): 15.0 mm
    elseif disktype == "E7R8"
        # set values for E7R8 ThinGap RRDEs
        # https://www.pineresearch.com/shop/products/electrodes/rrde/classic-ptfe-rrde/e7r8/
        # Collection Efficiency: 22%
        self.r_1 = 2.285 * mm
        self.r_2 = 2.465 * mm
        self.r_3 = 2.69 * mm
        self.r_4 = 5 * mm
    elseif disktype == "AlberyA"
        self.r_1 = 0.348 * cm
        self.r_2 = 0.386 * cm
        self.r_3 = 0.4375 * cm
        self.r_4 = 0.5 * cm
    elseif disktype == "AlberyB"
        self.r_1 = 0.3888 * cm
        self.r_2 = 0.3988 * cm
        self.r_3 = 0.4445 * cm
        self.r_4 = 0.5 * cm
    elseif disktype == "AlberyC"
        self.r_1 = 0.3672 * cm
        self.r_2 = 0.3763 * cm
        self.r_3 = 0.4369 * cm
        self.r_4 = 0.5 * cm
    elseif disktype == "AlberyD"
        self.r_1 = 0.3732 * cm
        self.r_2 = 0.3817 * cm
        self.r_3 = 0.492 * cm
        self.r_4 = 0.5 * cm
    elseif disktype == "AlberyE"
        self.r_1 = 0.25 * cm
        self.r_2 = 0.276 * cm
        self.r_3 = 0.359 * cm
        self.r_4 = 0.5 * cm
    elseif disktype == "AlberyF"
        self.r_1 = 0.25 * cm
        self.r_2 = 0.2738 * cm
        self.r_3 = 0.3592 * cm
        self.r_4 = 0.5 * cm
    else
        error("missing disktype")
    end
    #    self.w_gap=1.0e-2*self.r_4
    # We don't need this gap with real outflow BC.
    self.w_gap = 0.0
    self.h_cyl = 500 * μm
    return self
end
