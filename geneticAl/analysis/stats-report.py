#!/usr/bin/env python

import os
import sys

import SimpleUtils as utils
from Efficiency import *

efficient = Efficiency()

if __name__ == "__main__":

    monodatafile = open("monomer-zindo.txt", 'r')
    homos = {}
    lumos = {}
    for line in monodatafile:
        dataList = line.split(' ')
        if len(dataList) < 5:
            continue
        # (number, smiles, homo, lumo, oscstrs)
        smiles = dataList[1]
        homos[smiles] = float(dataList[2])
        lumos[smiles] = float(dataList[3])

    dimerdata = open("../dims_and_tets/alldimerDB.txt")
    dimHOMO = {}
    dimTrans = {}
    dimEff = {}
    dimHab = {}
    # skip line with column names
    dimerdata.readline()
    for line in dimerdata:
        dataList = line.split('"')
        if len(dataList) < 3:
            continue
        smiles = dataList[1]
        # smiles will include codes, like c(s1)c(S(=O)(=O)C=C2)c2c1~C#CC#C_2_1
        thisJson = dataList[5]
        (homo, trans) = utils.getHplusBG(thisJson)
        dimHOMO[smiles] = float(homo)
        dimTrans[smiles] = float(trans)
        dimEff[smiles] = efficient.zindoEff(float(homo), float(trans))
        mon1 = smiles[0:-4].split('~')[0]
        mon2 = smiles[0:-4].split('~')[1]
        print mon1, mon2
        # dimHab = 

    tetradata = open("../dims_and_tets/alltetramerDB.txt")
    tetradata.readline()
    for line in tetradata:
        dataList = line.split('"')
        smiles = dataList[1]
        thisJson = dataList[5]
        (homo, trans) = utils.getHplusBG(thisJson)
        print smiles, float(homo), float(trans), efficient.zindoEff(float(homo), float(trans)), dimHOMO[smiles], dimTrans[smiles], dimEff[smiles], dimHab[smiles]
