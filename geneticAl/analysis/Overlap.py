import os
import math
import sys

class Efficiency(object):
    def __init__(self):
        wavelen, spect = self.readdata()
        # Spect (W per m^2 per nm)
        # Wavelen, converting to m (from nm)
        wavelen = [x / 1.0E9 for x in wavelen]
        self.wavelen = wavelen
        self.spect = spect

        h = 6.626E-34 # J.s
        c = 3E8 # m per s

        spectrum_num = 1

        # phi = lambda * H / (hc)
        phi = [lamb * sp / (h * c) for lamb, sp in zip(self.wavelen,
                                                       self.spect[spectrum_num])]
        # Integrate
        area = [0]
        delta = [0]
        for j in range(1, len(phi)-1):
            delta_area = 0.5 * (phi[j] + phi[j+1]) * (self.wavelen[j+1] - self.wavelen[j])
            area.append(delta_area + area[-1])
            delta.append(delta_area)
        self.area = area
        delta.append(0)
        self.delta = delta

    def readdata(self, filename="AM0AM1_5_PVCDROM.txt"):
        relpath = os.sep.join(__file__.split(os.sep)[:-1])
        inputfile = open(os.path.join(relpath, filename), "r")
        header = inputfile.next()
        header = inputfile.next()
        wavelen = []
        spect = [[], [], []]
        maxLines = 0
        for line in inputfile:
            broken = map(float, line.split())
            wavelen.append(broken[0])
            for i in range(3):
                spect[i].append(broken[i+1])
            maxLines += 1

        self.maxLines = maxLines
        return wavelen, spect

    def gaussianFraction (self, x, center, width):
        # FWHM = 2(sqrt(2 ln 2)) * c
        c = width / (2.0* math.sqrt(2.0 * math.log(2.0)))
        return math.exp(-1.0*math.pow(x - center,2.0) / (2.0 * math.pow(c,2.0)))

    def efficiency(self, homo, bandgap, pcbm_lumo, verbose=False, cutoff=True):
        """HOMO,  bandgap in eV
        pcbm_lumo = -4.61 # Our calcs
        pcbm_lumo = -4.3 # Heeger
        """

        if bandgap == 0:
            # ZINDO calculation may give this value if it fails
            return 0
        if cutoff and homo + bandgap < (pcbm_lumo + 0.3):
            # => LUMO not high enough above PCBM to inject
            return 0

        Voc = abs(homo) - abs(pcbm_lumo) - 0.3
        if verbose:
            print "Voc", Voc

        q = 1.602E-19 # Coulumb (per electron)

        j = 0
        while (self.wavelen[j] * 1E9)< (1240/float(bandgap)):
            j += 1
        tot = self.area[j-1]
        if j==0:
            tot = 0

        tot *= 1.E9
        current = tot * q
        eff = current * 0.65 * Voc * 0.65
        eff /= 10.
        if verbose:
            print "# photons", tot, "Charge", current , "I (mA/cm2)", current * 0.1
            print "Efficiency", eff
        return eff

    def overlapCurrent(self, pBandgap, nBandgap, width):
        Voc = pBandgap - 0.3

        q = 1.602E-19 # Coulomb (per electron)

        # integrate the fraction of photons absorbed
        # across the entire spectrum
        total = 0.0
        for j in range(len(self.wavelen)):
            nm = self.wavelen[j] * 1.E9
            eV = 1240.0 / float(nm)
            pFraction = self.gaussianFraction(eV, pBandgap, width)
            nFraction = 0.0
            if (nBandgap >= pBandgap):
                nFraction = self.gaussianFraction(eV, nBandgap, width)
            fraction = pFraction + nFraction - pFraction*nFraction
            total += self.delta[j]*fraction
        total *= 1.E9
        current = total * q
        eff = current * 0.65 * Voc * 0.65
        eff /= 10.
        return current, eff

if __name__ == "__main__":
    eff = Efficiency()
    maxBG1 = 0.0
    maxBG2 = 0.0
    maxPower = 0.0
    width = 0.3
    for bg1 in range(10, 101):
        for bg2 in range(bg1, 101):
            pBg = float(bg1) / 25.0
            nBg = float(bg2) / 25.0
            if (nBg >= pBg):
                (current, power) = eff.overlapCurrent(pBg, nBg, width)
                print pBg, nBg, power, current
