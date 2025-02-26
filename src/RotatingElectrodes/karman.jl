"""
Data and methods for Karman RRDE velocity field
"""

mutable struct KarmanData
    omega
    nu

    A
    B
    a
    b
    alpha
    zeta_match

    Ap2
    Ap3
    Ap4
    Ap5
    Ap6

    Bp2
    Bp3
    Bp4
    Bp5
    Bp6

    ap2
    bp2
    ap4
    bp4

    alphap2
    alphap3
    alphap4
    alphap5
    alphap6
    alphap7
    alphap8
    alphap9
    alphap10
    alphap11

    zscale
    vzscale


    KarmanData(omega, nu) = KarmanData!(new(), omega, nu)
end

function KarmanData!(self, omega, nu)
    self.A = 0.919534683883937566757338594366182601566304873337718167234404388
    self.B = 1.185086462483026780888308224411023601905589737293560556689515672
    self.a = 0.5182057757142870936378111008096608891479003761353356941327685779
    self.b = -0.6300187813232994755200671106874782719915356089045425997266893724
    self.alpha = 0.8885668656643239309984114985119025125043523628112911617236627481
    self.zeta_match = 0.50272492317505845793590424364083446562290191650390625

    self.Ap2 = self.A * self.A
    self.Ap3 = self.A * self.A * self.A
    self.Ap4 = self.A * self.A * self.A * self.A
    self.Ap5 = self.A * self.A * self.A * self.A * self.A
    self.Ap6 = self.A * self.A * self.A * self.A * self.A * self.A

    self.Bp2 = self.B * self.B
    self.Bp3 = self.B * self.B * self.B
    self.Bp4 = self.B * self.B * self.B * self.B
    self.Bp5 = self.B * self.B * self.B * self.B * self.B
    self.Bp6 = self.B * self.B * self.B * self.B * self.B * self.B

    self.ap2 = self.a * self.a
    self.bp2 = self.b * self.b
    self.ap4 = self.a * self.a * self.a * self.a
    self.bp4 = self.b * self.b * self.b * self.b

    self.alphap2 = self.alpha * self.alpha
    self.alphap3 = self.alphap2 * self.alpha
    self.alphap4 = self.alphap3 * self.alpha
    self.alphap5 = self.alphap4 * self.alpha
    self.alphap6 = self.alphap5 * self.alpha
    self.alphap7 = self.alphap6 * self.alpha
    self.alphap8 = self.alphap7 * self.alpha
    self.alphap9 = self.alphap8 * self.alpha
    self.alphap10 = self.alphap9 * self.alpha
    self.alphap11 = self.alphap10 * self.alpha

    self.omega = omega
    self.nu = nu
    self.vzscale = sqrt(nu * omega)
    self.zscale = sqrt(omega / nu)

    #     println("KarmanData: omega=$(self.omega), nu=$(self.nu),zscale=$(self.zscale),vzscale=$(self.vzscale)")
    return self
end

function fgh(self, z0, linear = false)
    z = z0 * self.zscale
    f = 0.0
    g = 0.0
    h = 0.0
    if linear
        f = self.a * z
        g = 1.0 + self.b * z
        h = -self.a * z * z
    else
        if (z < self.zeta_match)
            f = (
                self.a - (
                    1.0 / 2.0 +
                        (self.b / 3.0 + (self.b * self.b / 12.0 + (self.a / 60.0 - 1.0 / 360.0 * (1.0 - 4.0 * self.a * self.b) * z) * z) * z) * z
                ) *
                    z
            ) * z
            g = 1.0 + (
                self.b + (
                    self.a / 3.0 + (
                        (1.0 / 12.0) * (-1.0 + self.a * self.b) -
                            (self.b / 15.0 - 1.0 / 90.0 * (-self.a * self.a - 2.0 * self.b * self.b) * z) * z
                    ) *
                        z
                ) *
                    z * z
            ) * z
            h = (
                -self.a + (
                    1.0 / 3.0 +
                        (self.b / 6.0 + (self.b * self.b / 30.0 + (self.a / 180.0 + ((-1.0 + 4.0 * self.a * self.b) * z) / 1260.0) * z) * z) * z
                ) *
                    z
            ) * z * z
        else
            ep1 = exp(-self.alpha * z)
            f = (
                (
                    (
                        (
                            (
                                -(1971.0 * self.Ap6 + (2825.0 * self.Ap4 + (889.0 * self.Ap2 + 35.0 * self.Bp2) * self.Bp2) * self.Bp2) / 86400.0 * ep1 / self.alphap2 -
                                    (-61.0 * self.Ap4 - (74.0 * self.Ap2 + 13.0 * self.Bp2) * self.Bp2) * self.A / 1152.0
                            ) *
                                ep1 / self.alphap2 -
                                (17.0 * self.Ap4 + (18.0 * self.Ap2 + self.Bp2) * self.Bp2) / 144.0
                        ) *
                            ep1 / self.alphap2 -
                            (-self.Ap3 - self.A * self.Bp2) / 4
                    ) *
                        ep1 / self.alphap2 -
                        (self.Ap2 + self.Bp2) / 2
                ) *
                    ep1 / self.alphap2 +
                    self.A
            ) *
                ep1
            g = (
                (
                    (
                        (
                            -(-65.0 * self.Ap4 - (82.0 * self.Ap2 + 17.0 * self.Bp2) * self.Bp2) * self.A * self.B / 5400.0 * ep1 / self.alphap2 -
                                (53.0 * self.Ap4 + (58.0 * self.Ap2 + 5.0 * self.Bp2) * self.Bp2) * self.B / 1920.0
                        ) *
                            ep1 / self.alphap2 -
                            (-self.Ap2 - self.Bp2) * self.A * self.B / 18.0
                    ) *
                        ep1 / self.alphap2 -
                        (self.Ap2 + self.Bp2) * self.B / 12.0
                ) *
                    ep1 * ep1 / self.alphap4 +
                    self.B
            ) *
                ep1
            h = (
                (
                    (
                        (
                            (
                                -(1971.0 * self.Ap6 + (2825.0 * self.Ap4 + (889.0 * self.Ap2 + 35.0 * self.Bp2) * self.Bp2) * self.Bp2) / 259200.0 * ep1 / self.alphap2 +
                                    self.A * (61.0 * self.Ap4 + (74.0 * self.Ap2 + 13.0 * self.Bp2) * self.Bp2) / 2880.0
                            ) *
                                ep1 / self.alphap2 -
                                (17.0 * self.Ap4 + (18.0 * self.Ap2 + self.Bp2) * self.Bp2) / 288.0
                        ) *
                            ep1 / self.alphap2 +
                            self.A * (self.A * self.A + self.B * self.B) / 6
                    ) *
                        ep1 / self.alphap2 -
                        (self.Ap2 + self.Bp2) / 2
                ) *
                    ep1 / self.alphap2 +
                    2 * self.A
            ) *
                ep1 / self.alpha -
                self.alpha
        end
    end
    return f, -g, h
end


function fKarman(self, r, z, linear = false)
    f, g, h = fgh(self, z, linear)
    v_r = self.omega * f * r
    v_z = self.vzscale * h
    return v_r, v_z
end
