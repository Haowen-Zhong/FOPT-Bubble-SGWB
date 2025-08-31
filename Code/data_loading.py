"""
该程序用来读取数据，并且存储为Delta_for_compare方便后续的区域计算。
要注意的是sigma=index+1
E.g.
	sigma(Delta_for_compare[0]=0+1)
"""

import numpy as np
i = np.arange(1, 101, 1)  # i 取值 1-100
Total_Delta = []#用来存储所有的数据
for index in i:
	data = np.genfromtxt(f"./Data/区域判定/sigma={index}.0.txt")
	Total_Delta.append(data)

Delta_for_compare = []
k_for_compare = []
for i in range(len(Total_Delta)):
	Delta_for_compare.append(Total_Delta[i][:,1])
	k_for_compare.append(Total_Delta[i][:,0])
