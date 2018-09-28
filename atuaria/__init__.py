import pandas as pd
import numpy as np


def tabuas():
    tabuas = pd.read_csv('https://raw.githubusercontent.com/lincolnfrias/dados/master/tabuas-de-vida.csv')
    return tabuas

def sv_vit(i, idade, b, qx):
    n = tabuas.idade.max() - idade
    px = 1 - qx.values
    serie = np.arange(1, n+1)
    v = 1/(i+1)**serie
    vp2 = (1/(i+1)**2)**serie
    qxx = qx[(idade):(idade+n)]
    pxx = np.cumprod(px[(idade):(idade+n-1)])
    pxx = np.insert(pxx, 0, 1)
    Ax = b * np.sum(v*pxx*qxx)
    Ax2 = b * np.sum(vp2*pxx*qxx)
    Var = (Ax2 - Ax**2)*b
    resultado = round(float(Ax), 1), round(float(Ax2), 1), round(float(Var), 1)
    return resultado

def sv_temp(i, idade, n, b, qx):
    px = 1 - qx.values
    serie = np.arange(1, n+1)
    v = 1/(i+1)**serie
    vp2 = (1/(i+1)**2)**serie
    qxx = qx[(idade):(idade+n)]
    pxx = np.cumprod(px[(idade):(idade+n-1)])
    pxx = np.insert(pxx, 0, 1)
    Ax = b * np.sum(v*pxx*qxx)
    Ax2 = b * np.sum(vp2*pxx*qxx)
    Var = (Ax2 - Ax**2)*b
    resultado = round(float(Ax), 1), round(float(Ax2), 1), round(float(Var), 1)
    return resultado
