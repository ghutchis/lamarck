import os
# Testing changes

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
        for j in range(1, len(phi)-1):
            delta_area = 0.5 * (phi[j] + phi[j+1]) * (self.wavelen[j+1] - self.wavelen[j])
            area.append(delta_area + area[-1])
        self.area = area

    def readdata(self, filename="AM0AM1_5_PVCDROM.txt"):
        relpath = os.sep.join(__file__.split(os.sep)[:-1])
        inputfile = open(os.path.join(relpath, filename), "r")
        header = inputfile.next()
        header = inputfile.next()
        wavelen = []
        spect = [[], [], []]
        for line in inputfile:
            broken = map(float, line.split())
            wavelen.append(broken[0])
            for i in range(3):
                spect[i].append(broken[i+1])

        return wavelen, spect

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

    def zindoEff(self, homo, bandgap):
        # HOMOs need no scaling, fit is 1.01 + ~0.4
        # PCBM LUMO = -6.7796 eV + (0.883 * 2.1886 + 0.530)
        #           = -4.317 eV
#        scaledGap = 0.883 * float(bandgap) + 0.530
#        return self.efficiency(homo, scaledGap, -4.3171)
#        return self.efficiency(homo, bandgap, -4.591)
        return self.efficiency(homo, bandgap, -4.61)
        # Numbers for PM6 -- currently unused
        # PCBM LUMO = -6.82998 eV + (1.1522*2.2187 + 0.0536)
        #           = -4.2736 eV
        # scaledGap = 1.1522*float(bandgap) + 0.0536
        # return self.efficiency(homo, scaledGap, -4.2736)

    def b3lypEff(self, homo, bandgap):
        # PCBM HOMO = 1.3023 * (-5.77368eV) + 0.481 = -7.0381 eV
        # Calibration from Zhan & Dixon, J. Phys. Chem. A (2003) 107 p. 4184.
        # PCBM LUMO = -7.0381 eV + (0.7397 * 2.5684eV + 0.8461)
        #           = -4.2922
        scaledHOMO = 1.3023 * float(homo) + 0.481
        scaledGap =  0.7397 * float(bandgap) + 0.8461
        return self.efficiency(scaledHOMO, scaledGap, -4.2922)

    def unittest(self):
        lumo = -3.92
        bandgap = 1.74
        homo = lumo - bandgap

        assert abs(self.efficiency(homo, bandgap, -4.3)-9.576) < 0.01, self.efficiency(homo, bandgap, -4.3)
        assert self.efficiency(-7.0 - 4.9, 4.9, -4.61) == 0, self.efficiency(-7.0 - 4.9, 4.9, -4.61)

if __name__ == "__main__":
    eff = Efficiency()
    eff.unittest()