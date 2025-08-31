#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Energy_density_spectrum import *
import sys

jobnumber = int(sys.argv[1])

data = open(f"./ratio_list_k_0.01.txt","w")
data.write("#      sigma     ratio\n")
sigma_list = np.arange(2,101,1)
ratio_list = np.arange(0.5,3.1,0.1)
for sigma in sigma_list:
    value = []
    for ix in range(len(ratio_list)):
        Delta_d = spectrum_d([0.01], sigma, ratio_list[ix])
        value.append(Delta_d)
        if len(value)>=2:
            if value[-1]<=value[-2]:
                data.write(f"{sigma}  {ratio_list[ix-1]}\n")
                break
data.close()