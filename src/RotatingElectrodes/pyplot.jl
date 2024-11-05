"""
   $(SIGNATURES)

    Plot CVs for a number of rotation frequencies.
    For each frequency, plot  disc curret and ring current vs. disk voltage.
"""
function pyplot_cv(freqs,vdisks,idisks,irings;
                   file="cv.pdf",
                   clear=true,
                   size_inches=(7,9)
                   )
    fig=PyPlot.figure(1)
    if clear
        PyPlot.clf()
    end
    fig.set_size_inches(size_inches...)
    PyPlot.grid(true,color="gray", linestyle="dotted")
    PyPlot.subplot(211)
    PyPlot.xlabel("Δϕ/V")
    PyPlot.ylabel("Disk current/mA")
    PyPlot.grid()
    nfreq=length(freqs)
    hfreq=1.0/(nfreq-1)
    xfreq=0.0
    for ifreq=1:length(freqs)
        freq=freqs[ifreq]
        idisk=idisks[ifreq]
        vdisk=vdisks[ifreq]
        lhalf=argmin(vdisk)
        PyPlot.plot(vdisk[1:lhalf],idisk[1:lhalf]./mA,label="f=$(freq) Hz",color=(xfreq,0,1.0-xfreq))
        PyPlot.plot(vdisk[lhalf:end],idisk[lhalf:end]./mA,"--",color=(xfreq,0,1.0-xfreq))
        PyPlot.legend(loc="upper right", fontsize="small", fancybox=false)
        xfreq+=hfreq
    end
    
    PyPlot.subplot(212)
    PyPlot.xlabel("Δϕ/V")
    PyPlot.ylabel("Ring current/mA")
    PyPlot.grid()
    nfreq=length(freqs)
    hfreq=1.0/(nfreq-1)
    xfreq=0.0
    for ifreq=1:length(freqs)
        freq=freqs[ifreq]
        vdisk=vdisks[ifreq]
        iring=irings[ifreq]
        lhalf=argmin(vdisk)
        PyPlot.plot(vdisk[1:lhalf],iring[1:lhalf]./mA,label="f=$(freq) Hz",color=(xfreq,0,1.0-xfreq))
        PyPlot.plot(vdisk[lhalf:end],iring[lhalf:end]./mA,"--",color=(xfreq,0,1.0-xfreq))
        
        PyPlot.legend(loc="upper right", fontsize="small", fancybox=false)
        PyPlot.tight_layout()
        xfreq+=hfreq
    end
    #PyPlot.tight_layout()
    PyPlot.savefig(file)
end
