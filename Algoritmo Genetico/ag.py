# -*- coding: utf-8 -*-
"""
Created on Thu May 18 12:15:04 2023

@author: THIAGO MORENO FERNANDES
"""


import math
import numpy as np
import array
import random
from piece_code_ag import get_bits, get_float
import matplotlib.pyplot as plt

def mutacao(s):
    pos = random.choice(range(len(s)))
    r = list(s)
    r[pos] = '1' if r[pos] == '0' else '0'

    return ''.join(r)


def ag(n=100, nger=500, nrun=5, Pc=0.5, pm=0.2, l=32):

    # Alocacao das variaveis
    pop_run = np.zeros((nrun, nger, n))
    fit_run = np.zeros((nrun, nger, n))

    for run in range(nrun):
        # Gera a populacao aleatoria
        pop = [random.uniform(0, math.pi) for _ in range(n)]
        pop = np.float32(pop)

        for ger in range(nger):

            # Avalia a funcao aptidao (fitness) da populacao
            fit = pop + np.abs(np.sin(32*pop))

            pop_run[run, ger, :] = pop
            fit_run[run, ger, :] = fit

            # Probabilidade de selecao
            Psel = []
            for f in fit:
                Psel.append(f/sum(fit))
                
            # Criacao de n filhos
            nova_pop = np.zeros(n)
            for nf in range(int(n/2)):
                # Seleciona os individuos - Selecao por roleta
                Cr = random.choices(pop, weights=Psel, k=2)

                # Transforma os indivíduos selecionados em bits
                Cr_bits = []
                for j in range(2):
                    Cr_bits.append(get_bits(Cr[j]))

                # Cruzamento
                # Distribuicao de probabilidade aleatoria uniformemente distribuida entre 0 e 1
                Cruz = random.random()

                # Define um ponto aleatoriamente no cromossomo - string de genes para cruzamento
                loci = random.randint(1, l-1)

                FilhoA = Cr_bits[0][:loci] + Cr_bits[1][loci:]
                FilhoB = Cr_bits[1][:loci] + Cr_bits[0][loci:]

                if Cruz >= Pc:  # Não ocorre o cruzamento
                    FilhoA = Cr_bits[0]
                    FilhoB = Cr_bits[1]

                # Mutacao do gene

                if random.random() <= pm:
                    FilhoA = [mutacao(s) for s in FilhoA]

                if random.random() <= pm:
                    FilhoB = [mutacao(s) for s in FilhoB]

                FilhoA = get_float(FilhoA)
                FilhoB = get_float(FilhoB)
                
                
                #Elimina os numeros NaNs trocando-os por pi
                if np.isnan(FilhoA):
                    FilhoA = np.pi
                if np.isnan(FilhoB):
                    FilhoB = np.pi

                nova_pop[nf*2:nf*2+2] = np.float32([FilhoA, FilhoB])

                #Se n for impar, um novo membro da populacao e descartado aleatoriamente
            if n % 2 > 0:
                nova_pop = np.delete(nova_pop, random.randint(n), 0)

            # Controle do dominio
            for e in range(n):
                if nova_pop[e] < 0:
                    nova_pop[e] = random.uniform(0, math.pi)
                elif nova_pop[e] > np.pi:
                    nova_pop[e] = random.uniform(0, math.pi)

            pop = nova_pop

    return pop_run, fit_run

# Executa o algoritmo genetico

n = 200
nger = 500
nrun = 5
Pc = 0.02
pm = 0.01
l = 32

pop_run, fit_run = ag(n=n, nger=nger, nrun=nrun, Pc=Pc, pm=pm, l=l)

"""
Plot da média da função fitness
"""
# Calcula a aptidao média de uma geracao para cada execucao
fit_media = np.mean(fit_run, axis=2, dtype=np.float64)

plt.rc('font', size=24)     
plt.rc('axes', titlesize=24) 
linestyle = [
    'solid',   # '-' 
    'dotted',  # ':'
    'dashed',  # '--'
    'dashdot', # '-.'
    ]
leg = []

for run in range(nrun):
    leg.append('Execução' + str(run))
    plt.plot(np.arange(nger), fit_media[run,:], linestyle=linestyle[run%len(linestyle)], linewidth=2.5, alpha=1-(run/nrun)*0.6)

# Eixos e legendas
plt.xlabel('Geracao')
plt.ylabel('Média da Aptidão')
plt.legend(leg)
plt.title(f'P_m = {pm}, P_c = {Pc}, População: {n}, Geração: {nger}')

# Salva a figura
fig = plt.gcf()
fig.set_size_inches(16, 9)
fig.savefig(f'{n}Execucao_{n}pop_{nger}gen_{Pc}Pc.png', dpi=300)

plt.show()

"""
Plota a população na função para a rodada 0
"""
geracao = [0, round(10-1), round(200-1), round(400-1), nger-1]
pop = pop_run[0, geracao, :]
fit = fit_run[0, geracao, :]

x = np.linspace(0, np.pi, num=1000)
y = x + np.abs(np.sin(32*x))
color = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ,'#1f77b4']
for idx, ger in enumerate(geracao):
    plt.figure(idx+2)
    plt.plot(x, y, linewidth=2.5)
    plt.scatter(pop[idx, :], fit[idx, :], c=color[idx])
    plt.legend(['Funcao', f'Ger {ger+1}'])
    plt.title(f'Geracao {ger+1}, P_m = {pm}, P_c = {Pc}, Populacao: {n}')
    fig = plt.gcf()
    fig.set_size_inches(16, 9)
    fig.savefig(f'{n}pop_ger{ger+1}_{Pc}Pc_aptidao.png', dpi=300)


plt.show()
