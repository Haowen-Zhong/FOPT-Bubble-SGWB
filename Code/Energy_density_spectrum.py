import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
from scipy.integrate import tplquad
from scipy.optimize import bisect
from scipy.integrate import quad


sin = np.sin
cos = np.cos
pi = np.pi
exp = np.exp
# 被积函数分块编码

def Prob(eta, sigma):
    return -4*pi/3/sigma**4*exp(eta-sigma)*(eta-sigma)**3


def False_Prob(duration,sigma,epsilon):
    return exp(quad(Prob,sigma,sigma+duration,args=(sigma))[0])-10**(-epsilon)

def tau(sigma, ratio):
    dur = bisect(False_Prob, a = 1, b = 30,args=(sigma,323,),maxiter=1000)
    dur = ratio*dur
    return dur

def Ixy(r, xi, Xi, sigma):
    """
    为了计算两时空点过去光锥内没有Bubble凝结的概率，我们需要计算出两者过去光锥所包含的
    四维共动体积，利用e^{-I(x,y)}即可得到概率P(x,y)
    Xi: (eta_x+eta_y)/2 这里也可以看出Xi_M=eta_f=sigma
    xi: eta_x-eta_y
    sigma: 用来衡量相变时标和宇宙膨胀时标。在我们的研究中希望sigma~O(1)
    """
    r2 = r*r
    r4 = r2*r2
    xi2 = xi*xi
    Xi2 = Xi*Xi
    PreFactor = exp(Xi-sigma)
    sigma2 = sigma*sigma
    sigma4 = sigma2*sigma2
    term1 = 96*PreFactor*exp(-0.5*xi)*r
    term2 = term1 * exp(xi)
    term3 = 24*PreFactor*exp(-0.5*r)*(-r*(4+r)+xi2)
    term4 = r4-4*r*(24*(1+Xi-sigma)+3*xi2*(1+Xi-sigma)+4*(Xi-sigma)*\
            (Xi-sigma)*(3+Xi-sigma))-12*xi2*(2+2*Xi+Xi2-2*(1+Xi)*sigma+\
            sigma2)-3*r2*(xi2+4*(2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2))
    LongTerm = term1 + term2 + term3 + term4
    Result = 1/(12*r*sigma4) * pi * LongTerm
    return Result


def F0(r, xi, Xi, sigma):
    """
    计算F0所用函数。
    """
    r2 = r*r
    r4 = r2*r2
    xi2 = xi*xi
    Xi2 = Xi*Xi
    sigma2 = sigma*sigma
    Result = 1/(16)*(r2-xi2)*(r2-xi2)*(32*exp(Xi-0.5*r-sigma)*(12+r*(6+r))-\
             (r4+384*(1+Xi-sigma)-8*r2*(2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2)+\
             16*(Xi-sigma)*(Xi-sigma)*(12+4*Xi+Xi2-2*(2+Xi)*sigma+sigma2)))
    return Result


def F1(r, xi, Xi, sigma):
    """
    计算F1所用代码
    """
    r2 = r*r
    r4 = r2*r2
    r6 = r4*r2
    xi2 = xi*xi
    Xi2 = Xi*Xi
    Xi3 = Xi2*Xi
    sigma2 = sigma*sigma
    Xi_sigma_2 = (Xi-sigma)*(Xi-sigma)
    Result = -1/(8)*(r2-xi2)*(16*exp(Xi-0.5*r-sigma)*(r2*(24+r*(12+r*(4+r)))-(120+r*(60+r*(12+r)))*\
             xi2)+(3*r6+1920*xi2+r4*(xi2-8*(2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2))+8*r2*(-48*(1+Xi-sigma)-\
             3*xi2*(2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2)-2*Xi_sigma_2*(12+4*Xi+Xi2-2*(2+Xi)*sigma+\
             sigma2))+80*xi2*(Xi-sigma)*(24+Xi3+Xi2*(4-3*sigma)-sigma*(12+(-4+sigma)*sigma)+Xi*(12+sigma*\
             (-8+3*sigma)))))
    return Result


def F2(r, xi, Xi, sigma):
    """
    计算F2所用代码
    """
    r2 = r*r
    r4 = r2*r2
    r6 = r2*r4
    r8 = r4*r4
    xi2 = xi*xi
    xi4 = xi2*xi2
    Xi2 = Xi*Xi
    Xi3 = Xi2*Xi
    Xi4 = Xi3*Xi
    sigma2 = sigma*sigma
    sigma3 = sigma2*sigma
    sigma4 = sigma3*sigma
    Xi_sigma_2 = (Xi-sigma)*(Xi-sigma)
    term1 = 96*(1+Xi-sigma)
    term2 = 2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2
    term3 = 4*Xi_sigma_2*(12+4*Xi+Xi2-2*(2+Xi)*sigma+sigma2)
    Result = -1/(16)*(3*r8+exp(Xi-0.5*r-sigma)*(-8*r4*(144+r*(72+r*(20+r*(4+r))))+16*r2*(720+r*\
             (360+r*(84+r*(12+r))))*xi2-8*(1680+r*(840+r*(180+r*(20+r))))*xi4)+560*xi4*(24+24*Xi+12*Xi2+\
             4*Xi3+Xi4-4*(6+Xi*(6+Xi*(3+Xi)))*sigma+6*(2+Xi*(2+Xi))*sigma2-4*(1+Xi)*sigma3+sigma4)+\
             2*r6*(xi2+4*term2)-120*r2*xi2*(term1+xi2*term2+term3)+3*r4*(xi4+4*term1+16*xi2*term2+4*term3))
    return Result


def j0(x):
    """
    0阶球贝塞尔函数j0
    为了防止出问题，认为设定x小于0.001时取极限值1n
    """
    if x < 0.001:
        return 1
    else:
        return sin(x)/x


def j1(x):
    """
    1阶球贝塞尔函数j1/x
    为了防止出问题，认为设定x小于0.001时取极限值1/3
    """
    if x < 0.001:
        return 1/3
    else:
        return (sin(x)-x*cos(x))/(x*x*x)


def j2(x):
    """
    2阶球贝塞尔函数j2/x^2
    为了防止出问题，认为设定x小于0.001时取极限值1/15
    """
    if x < 0.001:
        return 1/15
    else:
        return ((3-x*x)*sin(x)-3*x*cos(x))/(x*x*x*x*x)

def Bx(r,xi,Xi,sigma):
    r2 = r*r
    r4 = r2*r2
    xi2 = xi*xi
    Xi2 = Xi*Xi
    sigma2 = sigma*sigma
    return (r2-xi2)*(4*exp(Xi-0.5*r-sigma)*(12*xi+r*(6*xi+r*(2+r+xi)))+(r4-48*xi-8*xi*\
           (Xi-sigma)*(6+Xi2+Xi*(3-2*sigma)+(-3+sigma)*sigma)+2*r2*(xi*(1+Xi-sigma)-2*\
           (2+2*Xi+Xi2-2*(1+Xi)*sigma+sigma2))))/24
    
def By(r, xi, Xi, sigma):
    r2 = r*r
    r4 = r2*r2
    xi2 = xi*xi
    Xi2 = Xi*Xi
    sigma2 = sigma*sigma
    return (r2-xi2)*(4*exp(Xi-0.5*r-sigma)*(r2*(2+r)-(12+r*(6+r))*xi)+(r4+48*xi+8*xi*(Xi-sigma)*(6+Xi2+\
           Xi*(3-2*sigma)+(-3+sigma)*sigma)-2*r2*(4+xi+4*Xi+xi*Xi+2*Xi2-(4+xi+4*Xi)*sigma+2*sigma2)))/24
    
def Integrand_s(r,xi,Xi,k,sigma):
    """
    三重积分被积函数
    k: 共动波矢 co-moving wave vector
    r: x,y两时空点共动坐标差 r:=|x-y|
    """
    I = Ixy(r, xi, Xi, sigma)
    x = k*r
    LongTerm = j0(x)*F0(r, xi, Xi, sigma) + j1(x)*F1(r, xi, Xi, sigma) + j2(x)*F2(r, xi, Xi, sigma)
    Result = (Xi*Xi-xi*xi/4)*(Xi*Xi-xi*xi/4)*(Xi*Xi-xi*xi/4)/(r*r*r)*cos(k*xi)*exp(-I)*LongTerm
    return Result


def Integrand_d(r,xi,Xi,k,sigma):
    I = Ixy(r, xi, Xi, sigma)
    P = exp(-I)
    x = k*r
    Result = (Xi*Xi-xi*xi/4)*(Xi*Xi-xi*xi/4)*(Xi*Xi-xi*xi/4)*cos(k*xi)*P/(r*r*r*r)*Bx(r,xi,Xi,sigma)*By(r,xi,Xi,sigma)*j2(x)
    return Result


def spectrum_s(k, sigma, ratio):
    Delta_s = []
    time_i = time.time()
    duration = tau(sigma, ratio)
    for wavevec in k:
        result = tplquad(Integrand_s,  # 被积函数
                         sigma,  # Xi下限
                         sigma+duration,  # Xi上限
                         lambda Xi: 0,  # xi下限
                         lambda Xi: duration,  # xi上限
                         lambda Xi, xi: xi,  # r下限
                         lambda Xi, xi: 2*(Xi-sigma),  # r上限
                         args=(wavevec, sigma),
                         epsrel=0.005
                         )[0]
        result = 2/3*wavevec*wavevec*wavevec*result/(sigma*sigma*sigma*sigma)/(sigma+duration)**8
        Delta_s.append(result)
    Delta_s = np.array(Delta_s)
    time_f = time.time()
    print(f"计算结束，共花费{(time_f - time_i) / 60}分钟。")
    return Delta_s


def spectrum_d(k, sigma, ratio):
    Delta_d = []
    time_i = time.time()
    duration = tau(sigma, ratio)
    for wavevec in k:
        result = tplquad(Integrand_d,  # 被积函数
                         sigma,  # Xi下限
                         sigma+duration,  # Xi上限
                         lambda Xi: 0,  # xi下限
                         lambda Xi: duration,  # xi上限
                         lambda Xi, xi: xi,  # r下限
                         lambda Xi, xi: 2*(Xi-sigma),  # r上限
                         args=(wavevec, sigma),
                         epsrel=0.005
                         )[0]
        result = 24*pi*wavevec*wavevec*wavevec*result/(sigma*(sigma+duration))**8
        Delta_d.append(result)
    Delta_d = np.array(Delta_d)
    time_f = time.time()
    print(f"计算结束，共花费{(time_f - time_i) / 60}分钟。")
    return Delta_d
