# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:05:35 2023
@author: Thiago Moreno Fernandes
"""

print('Autor: Thiago Moreno Fernandes')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#inputs
i = 5000   # Numero de iterações
alpha = 0.1   #Taxa de aprendizado

theta_hat = np.zeros((i, 2))
custo_hat = np.zeros((i, 1))

theta_inicial = np.array((0, 0.1)) #theta inicial

def f_true (x):
    
    return 2 + 0.8 * x
    
# conjunto de dados {(x,y)}
xs = np. linspace (-3, 3, 100)
ys = np. array ( [ f_true (x) + np. random . randn () *0.5 for x in xs] )
m = len(xs)

#hipotese
def h(x, theta):
    
    htheta = theta[0]+ theta[1]*x
    
    return htheta


#funcao custo
def J(theta, xs , ys, m):
    
    custo = (1/(2*m))*np.sum((h(xs,theta)-ys)**2)
    
    return custo


#derivada parcial com respeito a theta [i]
def gradient(i, alpha, theta, xs, ys, m):
    
    for epocas in range(i):

            dtheta0 = np.array(h(xs,theta)-ys)
            dtheta1 = np.array((h(xs,theta)-ys)*xs)

            theta_hat0 = theta[0] - (alpha/m)*np.sum(dtheta0)
            theta_hat1 = theta[1] - (alpha/m)*np.sum(dtheta1)
            
            theta = np.array((theta_hat0, theta_hat1))
            
            theta_hat[epocas,:] = theta
            custo_hat[epocas,:] = J(theta, xs , ys, m)
                                             
    return theta,theta_hat,custo_hat


                                     
#plota no mesmo grafico : - o modelo / hipotese (reta) a reta original ( true function ) e os dados com ruido (xs , ys)"
def print_modelo (ax, theta, xs, ys, c='r', label=None):
    
    x = xs
    y = f_true(xs) #Função verdadeira
    y_hat = h(xs, theta)  #Função estimada pela regressão
    
    
    ax.scatter(x, ys, linewidths=1.0, c='gray', marker='.', label='Dados')
    ax.plot(x, y, linewidth=2.0, c='k', label='Função verdadeira')
    ax.plot(x, y_hat, linewidth=2.0, linestyle='--', c=c, label='Função estimada')
    ax.legend()
    
    fig.suptitle('Função verdadeira x Função estimada')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.rc('font', size=10)

theta, theta_hat, custo_hat = gradient(i=i, alpha=alpha, theta=theta_inicial, xs=xs, ys=ys, m=m)
fig, ax = plt.subplots()
print_modelo (ax, theta, xs, ys, c='r', label='label')

# Plota a função custo ao longo das iterações (épocas)

def plot_epocas (i, custo_hat, c='r', label=None):
    
    
    
    
    
        








