import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
sns.set_style("ticks")
DECIGO = np.genfromtxt("./Data/区域判定/Sensitivity_curve_DECIGO.txt")
BBO = np.genfromtxt("./Data/区域判定/Sensitivity_curve_BBO.txt")
LISA = np.genfromtxt("./Data/区域判定/Sensitivity_curve_LISA.txt")
TianQin = np.genfromtxt("./Data/区域判定/Sensitivity_curve_TianQin.txt")
B_DECIGO = np.genfromtxt("./Data/区域判定/Sensitivity_curve_B_DECIGO.txt")
EPTA = np.genfromtxt("./Data/区域判定/Sensitivity_curve_EPTA.txt")
NANOGrav = np.genfromtxt("./Data/区域判定/Sensitivity_curve_NANOGrav.txt")
SKA5 = np.genfromtxt("./Data/区域判定/Sensitivity_curve_SKA5.txt")
SKA10 = np.genfromtxt("./Data/区域判定/Sensitivity_curve_SKA10.txt")
SKA20 = np.genfromtxt("./Data/区域判定/Sensitivity_curve_SKA20.txt")
k_LISA = []
sensitivity_LISA = []
k_TianQin = []
sensitivity_TianQin = []
k_DECIGO = []
sensitivity_DECIGO = []
k_BBO = []
sensitivity_BBO = []
k_B_DECIGO = []
sensitivity_B_DECIGO = []
k_EPTA = []
sensitivity_EPTA = []
k_NANOGrav = []
sensitivity_NANOGrav = []
k_SKA5 = []
sensitivity_SKA5 = []
k_SKA10 = []
sensitivity_SKA10 = []
k_SKA20 = []
sensitivity_SKA20 = []
for item in DECIGO:
    k_DECIGO.append(item[0])
    sensitivity_DECIGO.append(item[1])
for item in BBO:
    k_BBO.append(item[0])
    sensitivity_BBO.append(item[1])
for item in B_DECIGO:
    k_B_DECIGO.append(item[0])
    sensitivity_B_DECIGO.append(item[1])
for item in LISA:
    k_LISA.append(item[0])
    sensitivity_LISA.append(item[1])
for item in TianQin:
    k_TianQin.append(item[0])
    sensitivity_TianQin.append(item[1])
for item in EPTA:
    k_EPTA.append(item[0])
    sensitivity_EPTA.append(item[1])
for item in NANOGrav:
    k_NANOGrav.append(item[0])
    sensitivity_NANOGrav.append(item[1])
for item in SKA5:
    k_SKA5.append(item[0])
    sensitivity_SKA5.append(item[1])
for item in SKA10:
    k_SKA10.append(item[0])
    sensitivity_SKA10.append(item[1])
for item in SKA20:
    k_SKA20.append(item[0])
    sensitivity_SKA20.append(item[1])
def senstivity_curves():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.loglog(k_DECIGO,sensitivity_DECIGO,"--",color="forestgreen",linewidth=1)
    plt.fill_between(k_DECIGO,sensitivity_DECIGO,1,facecolor="forestgreen",alpha=0.2)
    plt.loglog(k_BBO,sensitivity_BBO,"--",color="purple",linewidth=1)
    plt.fill_between(k_BBO,sensitivity_BBO,sensitivity_DECIGO,facecolor="purple",alpha=0.2)
    plt.loglog(k_B_DECIGO,sensitivity_B_DECIGO,"--",color="royalblue",linewidth=1)
    plt.fill_between(k_B_DECIGO,sensitivity_B_DECIGO,1,facecolor="royalblue",alpha=0.2)
    plt.loglog(k_LISA,sensitivity_LISA,"--",color="olivedrab",linewidth=1)
    plt.fill_between(k_LISA,sensitivity_LISA,1,facecolor="olivedrab",alpha=0.2)
    plt.loglog(k_TianQin,sensitivity_TianQin,"--",color="green",linewidth=1)
    plt.fill_between(k_TianQin,sensitivity_TianQin,1,facecolor="green",alpha=0.2)
    plt.loglog(k_EPTA,sensitivity_EPTA,"--",color="red",linewidth=1)
    plt.fill_between(k_EPTA,sensitivity_EPTA,1,facecolor="red",alpha=0.2)
    plt.loglog(k_NANOGrav,sensitivity_NANOGrav,"--",color="chocolate",linewidth=1)
    plt.fill_between(k_NANOGrav,sensitivity_NANOGrav,1,facecolor="chocolate",alpha=0.2)
    plt.loglog(k_SKA5,sensitivity_SKA5,"--",color="orange",linewidth=1)
    plt.fill_between(k_SKA5,sensitivity_SKA5,1,facecolor="orange",alpha=0.2)
    plt.loglog(k_SKA10,sensitivity_SKA10,"--",color="orange",linewidth=1)
    plt.fill_between(k_SKA10,sensitivity_SKA10,1,facecolor="orange",alpha=0.2)
    plt.loglog(k_SKA20,sensitivity_SKA20,"--",color="orange",linewidth=1)
    plt.fill_between(k_SKA20,sensitivity_SKA20,1,facecolor="orange",alpha=0.2)
    plt.text(6e-1,1e-11,"B-DECIGO",rotation=50,color="royalblue",fontsize=12)
    plt.text(1.2e-2,0.7e-14,"DECIGO",color="forestgreen",fontsize=12)
    plt.text(5,1e-14,"BBO",rotation=50,color="purple",fontsize=12)
    plt.text(3e-8,7e-9,"NG",rotation=50,color="chocolate",fontsize=12)
    plt.text(5e-9,1e-8,"EPTA",rotation=50,color="red",fontsize=12)
    plt.text(2e-5,2e-11,"TianQin",rotation=-40,color="green",fontsize=12)
    plt.text(9e-6,9e-13,"LISA",rotation=-40,color="olivedrab",fontsize=12)
    plt.text(1.2e-7,2e-11,"SKA",rotation=50,color="orange",fontsize=12)
    plt.text(5e-8,1e-13,"5 years",rotation=-35,color="orange",fontsize=8)
    plt.text(1.5e-8,1e-15,"10 years",rotation=-35,color="orange",fontsize=8)
    plt.text(3e-9,1e-17,"20 years",rotation=-35,color="orange",fontsize=8)
#plt.savefig("../Figure/PLI.pdf",dpi=500,bbox_inches="tight")
f_min = k_BBO[np.argmin(sensitivity_BBO)]
senstivity_min = np.min(sensitivity_BBO)