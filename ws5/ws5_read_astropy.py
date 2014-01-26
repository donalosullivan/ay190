#!/usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as pl

data = ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat")  
logsigma = np.array(np.log10(data["sigma*"]))
logM = np.array(data["logM"])

