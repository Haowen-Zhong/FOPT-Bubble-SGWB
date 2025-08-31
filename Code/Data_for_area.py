#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Energy_density_spectrum import *
import multiprocessing as mp
import sys

sigma_min = float(sys.argv[1])
sigma_max = float(sys.argv[2])


num_cores = int(mp.cpu_count())
print("本地计算机有: " + str(num_cores) + " 核心")
pool = mp.Pool(num_cores)

k = np.logspace(0, 1, 15)
sigma = np.arange(sigma_min,sigma_max,1)
ratio_list = np.genfromtxt("./ratio_list.txt")
for ix in range(len(sigma)):
	print(f"已经开始计算sigma={sigma[ix]}的情形")
	time_i = time.time()
	data = open(f"../Data/区域判定/sigma={sigma[ix]}.txt","w")
	data.write(f"#   sigma={sigma[ix]}\n")
	data.write("# k   Delta_s   Delta_d\n")

	results = [pool.apply_async(spectrum_s, args=(k, sigma[ix], ratio_list[int(sigma[ix]-2),1])), pool.apply_async(spectrum_d, args=(k, sigma[ix], ratio_list[int(sigma[ix]-2),1]))]
	Delta_s, Delta_d = [p.get() for p in results]
	for i in range(len(k)):
		data.write(f"{k[i]}  {Delta_s[i]+Delta_d[i]}\n")
	data.close()
	time_f = time.time()
	duration = time_f - time_i
	print(f"sigma={sigma[ix]}的情况已经计算完成，该次计算用时{duration/60}分钟")

