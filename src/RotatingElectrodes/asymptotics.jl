"""
Classical asymptotic and analytical expressios for limiting currents and collection efficiencies
"""

function idisk_levich(self,omega,nu,D)
    FaradayConstant=ph"AvogadroConstant"*ph"ElementaryCharge"
    Sc=nu/D
    return FaradayConstant*self.r_1*self.r_1*π*(omega*nu)^(1.0/2)*0.6205*Sc^(-2.0/3.)
end

function iring_levich(self,omega,nu,D)
    FaradayConstant=ph"AvogadroConstant"*ph"ElementaryCharge"
    return  (FaradayConstant*π*(1.0/(1.61))
             *D^(2.0/3.0)
             *nu^(-1.0/6.0)
             *omega^(1.0/2.0)
             *(self.r_3^3.0-self.r_2^3.0)^(2.0/3.0))
end

function idisk_newman(self,omega,nu,D)
    FaradayConstant=ph"AvogadroConstant"*ph"ElementaryCharge"
    Sc=nu/D
    return FaradayConstant*self.r_1*self.r_1*π*(omega*nu)^(1.0/2)*0.62048*Sc^(-2.0/3.)/(1.0+0.2980*Sc^(-1.0/3.)+0.14514*Sc^(-2.0/3.))
end

dlayer(nu,diff,omega)=1.61*(nu)^(1/6)*(diff)^(1/3)/sqrt(omega)

function coleff_albery(self)
    function F(theta)
        thetapd = theta^(1.0/3.0)
        return ( (sqrt(3.0)/(4.0*π))*log((1.0+thetapd)^3/(1.0+theta))
                 + (3.0/(2.0*π))*atan((2.0*thetapd-1.0)/sqrt(3.0))
                 + 1.0/4.0)
    end
    alpha = (self.r_2/self.r_1)^3 - 1.0
    beta = (self.r_3*self.r_3*self.r_3 - self.r_2*self.r_2*self.r_2)/(self.r_1*self.r_1*self.r_1)
    return (1.0
            -F(alpha/beta)
            +beta^(2.0/3.0)*(1.0-F(alpha))
            -(1+alpha+beta)^(2.0/3.0)*(1-F((alpha/beta)*(1.0+alpha+beta))))
end
