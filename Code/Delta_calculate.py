#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Energy_density_spectrum import *
import sys

jobnumber = int(sys.argv[1])
sigma = float(sys.argv[2])
ratio = float(sys.argv[3])

data = open(f"./Energy_density/data_sigma={sigma}_ratio={ratio}.txt","w")
data.write(f"#   sigma={sigma}  ratio={ratio}\n")
data.write("# k Delta_tot\n")
plt.rc('text', usetex=True)
k = np.logspace(-2,1,30)
Delta_s = spectrum_s(k, sigma, ratio)
Delta_d = spectrum_d(k, sigma, ratio)
for i in range(len(k)):
	data.write(f"{k[i]}  {Delta_s[i]+Delta_d[i]}\n")
data.close()