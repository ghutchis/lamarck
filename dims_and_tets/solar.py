import pdb
import os
import numpy
import pylab
import matplotlib
import matplotlib.image as image

import sys
sys.path.append(os.path.join("..", "geneticAl"))

import Efficiency as Eff

##def readdata(filename="AM0AM1_5_PVCDROM.txt"):
##    inputfile = open(filename, "r")
##    header = inputfile.next()
##    header = inputfile.next()
##    wavelen = []
##    spect = [[], [], []]
##    for line in inputfile:
##        broken = map(float, line.split())
##        wavelen.append(broken[0])
##        for i in range(3):
##            spect[i].append(broken[i+1])
##        
##    return wavelen, spect

##def unittest():
##    
##    wavelen, spect = readdata()
##    # Spect (W per m^2 per nm)
##    # Wavelen, converting to m (from nm)
##    wavelen = [x / 1E9 for x in wavelen]
##
##    q = 1.602E-19 # Coulumb (per electron)
##    h = 6.626E-34 # J.s
##    c = 3E8 # m per s
##
##    for i in range(3):
##        # phi = lambda * H / (hc)
##        phi = [lamb * sp / (h * c) for lamb, sp in zip(wavelen, spect[i])]
##
##        # Integrate
##        tot = 0
##        silicon = 0
##        for j in range(len(phi)-1):
##            area = 0.5 * (phi[j] + phi[j+1]) * (wavelen[j+1] - wavelen[j])
##            tot += area
##            if (wavelen[j] * 1E9)< 1107: # For all wavelengths below the bandgap of Si
##                silicon += area
##
##        tot *= 1E9
##        silicon *= 1E9
##        print "# photons", tot, "Charge", tot * q, "I (mA/cm2)", tot*q*0.1
##        print "Si photons", silicon, "Charge", silicon * q, "I (mA/cm2)", silicon*q*0.1
##
##        if i==1:
##            assert abs(tot*q*0.1 - 69) < 1 # 69mA/cm2
##            assert abs(silicon*q*0.1 - 43.8) < 0.1 # 43.8A/cm2
##    
##
##class Efficiency(object):
##    def __init__(self):
##        wavelen, spect = readdata()
##        # Spect (W per m^2 per nm)
##        # Wavelen, converting to m (from nm)
##        wavelen = [x / 1.0E9 for x in wavelen]
##        self.wavelen = wavelen
##        self.spect = spect
##        
##    def efficiency(self, lumo, bandgap, verbose=True):
##        # lumo, bg in eV
##
##        homo = lumo - bandgap
##        Voc = abs(homo) - 4.3 - 0.3
##        if verbose:
##            print "Voc", Voc
##
##        q = 1.602E-19 # Coulumb (per electron)
##        h = 6.626E-34 # J.s
##        c = 3E8 # m per s
##
##        spectrum_num = 1
##
##        # phi = lambda * H / (hc)
##        phi = [lamb * sp / (h * c) for lamb, sp in zip(self.wavelen,
##                                                       self.spect[spectrum_num])]
##
##        # Integrate
##        tot = 0
##        for j in range(len(phi)-1):
##            area = 0.5 * (phi[j] + phi[j+1]) * (self.wavelen[j+1] - self.wavelen[j])
##            if (self.wavelen[j] * 1E9)< (1240/float(bandgap)): # For all wavelengths below the bandgap of Si
##                tot += area
##        tot *= 1E9
##        current = tot * q
##        eff = current * 0.65 * Voc * 0.65
##        eff /= 10.
##        if verbose:
##            print "# photons", tot, "Charge", current , "I (mA/cm2)", current * 0.1
##            print "Efficiency", eff
##        return eff
##
##class Efficiency_b(object):
##    def __init__(self):
##        wavelen, spect = readdata()
##        # Spect (W per m^2 per nm)
##        # Wavelen, converting to m (from nm)
##        wavelen = [x / 1.0E9 for x in wavelen]
##        self.wavelen = wavelen
##        self.spect = spect
##
##        h = 6.626E-34 # J.s
##        c = 3E8 # m per s
##
##        spectrum_num = 1
##
##        # phi = lambda * H / (hc)
##        phi = [lamb * sp / (h * c) for lamb, sp in zip(self.wavelen,
##                                                       self.spect[spectrum_num])]
##        # Integrate
##        area = [0]
##        for j in range(1, len(phi)-1):
##            delta_area = 0.5 * (phi[j] + phi[j+1]) * (self.wavelen[j+1] - self.wavelen[j])
##            area.append(delta_area + area[-1])
##        self.area = area
##        
##    def efficiency(self, lumo, bandgap, verbose=True, cutoff=True, pcbm_lumo=-4.61):
####        pcbm_lumo = -4.61 # Our calcs
####        pcbm_lumo = -4.3 # Heeger
##        
##        if cutoff is True:
##            # lumo, bg in eV
##            if lumo < (pcbm_lumo + 0.3): # This is the lowest E PCBM excitation
##                return 0
##
##        homo = lumo - bandgap
##        
##        Voc = abs(homo) - abs(pcbm_lumo) - 0.3
##        if verbose:
##            print "Voc", Voc
##
##        q = 1.602E-19 # Coulumb (per electron)
##
##        j = 0
##        while j<len(self.wavelen)-1 and (
##              self.wavelen[j] * 1E9)< (1240/float(bandgap)):
##            j += 1
##        tot = self.area[j-1]
##        if j == 0:
##            tot = 0
##        tot *= 1.E9
##        current = tot * q
##        eff = current * 0.65 * Voc * 0.65
##        eff /= 10.
##        if verbose:
##            print "# photons", tot, "Charge", current , "I (mA/cm2)", current * 0.1
##            print "Efficiency", eff
##        return eff
##
##def test():
##    for effcalc in Efficiency(), Efficiency_b():
##        assert abs(effcalc.efficiency(-3.92, 1.74, verbose=False)-9.576) < 0.01
##        assert effcalc.efficiency(-7.0, 4.9, verbose=False) == 0
##
##def heeger():
##    effcalc = Efficiency_b()
##
##    im = image.imread(os.path.join("..", "Geoff-8Sep", "Efficiency.png"))
####    im[:, :, -1] = 0.0 # set the alpha channel
##    pylab.imshow(im, extent=(3.1, 1.0, -3.0, -4.0))
##    
##    alleff = []
##    delta = 0.05
####    delta = 0.5
##    lims = [[-2.0, -5.0], [5.0, 0.5]]
##    ranges = [numpy.arange(lims[0][0], lims[0][1] - delta, -delta),
##              numpy.arange(lims[1][0], lims[1][1]-delta, -delta)]
##    for lumo in ranges[0]:
##        eff = []
##        for bandgap in ranges[1]:
##            eff.append(effcalc.efficiency(lumo, bandgap, -4.3,
##                                          verbose=False, cutoff=False))
##        alleff.append(eff)
##    cs = pylab.contour(ranges[1], ranges[0], numpy.array(alleff), numpy.arange(1.0, 40, 1.0))
##
##    pylab.xlim(lims[1][0], lims[1][1])
##    pylab.ylim(lims[0][0], lims[0][1])
##    pylab.clabel(cs, inline=0, fontsize=10)
##    pylab.show()

def plot(pcbm_lumo, colorbar=True):
    """Plot the landscape of efficiency values

    This is imported by other scripts"""
    
    effcalc = Eff.Efficiency()

    alleff = []
    delta = 0.05

    lims = [[-6.8, -4.7], [2.5, 0.5]]
    lims = [[-11.5, -4.7], [5.3, 0.5]]
    ranges = [numpy.arange(lims[0][0], lims[0][1] + delta, delta),
              numpy.arange(lims[1][0], lims[1][1]-delta, -delta)]

    cutoff = pcbm_lumo + 0.3
    contours = [0, 0.0001, 0.01, 0.04, 0.1, 0.2, 0.4, 0.7, 1,2,3,4,5, 6, 8, 10, 12, 14]
    labels = ["0", "1E-4", "0.01", "0.04", "0.1", "0.2", "0.4", "0.7"]
    labels += ["%d" % x for x in contours if x >= 1.0]
    for homo in ranges[0]:
        eff = []
        for bandgap in ranges[1]:
            eff.append(effcalc.efficiency(homo, bandgap, pcbm_lumo,
                                          verbose=False, cutoff=False))
        alleff.append(eff)
    cs = pylab.contour(ranges[1], ranges[0], numpy.array(alleff), contours)
    pylab.plot([0, 6.5], [cutoff, -6.5 + cutoff], color="black")
    pylab.xlim(lims[1][1], lims[1][0])
    pylab.ylim(lims[0][0], lims[0][1])
    if colorbar:
        CB = pylab.colorbar(cs, shrink=0.8, extend='both', ticks=contours,
                            format=matplotlib.ticker.FixedFormatter(labels),
                            )
        CB.set_label("% Efficiency")
        return CB

def plot_heeger(pcbm_lumo, colorbar=True):
    """Plot the landscape of efficiency values

    This is imported by other scripts"""

    cutoff = pcbm_lumo + 0.3    
    effcalc = Eff.Efficiency()

    alleff = []
    delta = 0.05

    lims = [[-6.8, -4.7], [2.5, 0.5]]
    lims = [[cutoff+1.0, cutoff], [1.0, 3.1]]
    ranges = [numpy.arange(lims[0][0], lims[0][1] - delta, -delta),
              numpy.arange(lims[1][0], lims[1][1] + delta, delta)]
    print

    contours = [0, 0.5, 1,2,3,4,5, 6, 7, 8, 9, 10, 11]
    labels = ["0", "0.5"]
    labels += ["%d" % x for x in contours if x >= 1.0]
    for homo_plus_bg in ranges[0]:
        eff = []
        for bandgap in ranges[1]:
            homo = homo_plus_bg - bandgap
            eff.append(effcalc.efficiency(homo, bandgap, pcbm_lumo,
                                          verbose=False, cutoff=False))
        alleff.append(eff)
    cs = pylab.contour(ranges[1], ranges[0], numpy.array(alleff), contours)
##    pylab.plot([0, 6.5], [cutoff, -6.5 + cutoff], color="black")
    pylab.xlim(lims[1][1], lims[1][0])
    pylab.ylim(lims[0][0], lims[0][1])
    if colorbar:
        CB = pylab.colorbar(cs, shrink=0.8, extend='both', ticks=contours,
                            format=matplotlib.ticker.FixedFormatter(labels),
                            )
        CB.set_label("% Efficiency")
        return CB

def find_max(pcbm_lumo):
    """Find the maximum efficiency value for a particular
    pcbm_lumo"""
    
    effcalc = Eff.Efficiency()
    delta = 0.005
##    delta = 0.5
    lims = [[-8.5, -4.91], [2.5, 0.5]]
    ranges = [numpy.arange(lims[0][0], lims[0][1] + delta, delta),
              numpy.arange(lims[1][0], lims[1][1]-delta, -delta)]

    cutoff = pcbm_lumo + 0.3
    eff = []
    maxe = (0, 0, 0, 0)
    for homo in ranges[0]:
        lumo = cutoff
        bandgap = lumo - homo
        e = effcalc.efficiency(lumo, bandgap, False, False,
                                      pcbm_lumo=pcbm_lumo)
        eff.append(e)
        if e > maxe[0]:
            maxe = (e, homo, lumo, bandgap)
    print "Max", maxe

    pylab.plot(ranges[0], eff)
##    pylab.xlim(lims[1][1], lims[1][0])
##    CB = pylab.colorbar(cs, shrink=0.8, extend='both', ticks=contours)    

if __name__=="__main__":
##    find_max(-4.3)
##    find_max(-4.61)
##    pylab.show()

    plot_heeger(-4.61)
    pylab.show()
            
            
    
    
