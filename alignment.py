# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 09:16:41 2020

@author: 1
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


source_data = {'L': 150,
               'n': 999,
               'n+1': 1000,
               'm_B': 0.5,
               'ro': 206265
        }
source_data['ΔD'] = source_data['L'] / source_data['n+1']
#%% Полный створ


def full_target(source_data):
    m_i = []
    for i in range(source_data['n+1']):
        if i != 0:
            a = ((np.sqrt(2) * source_data['m_B'] * source_data['L'] * 10 ** 3 * 
                  i * (source_data['n+1'] - i)) /
                 (source_data['ro'] * source_data['n+1'] * np.sqrt(i ** 2 + (
                         source_data['n+1'] - i) ** 2)))
            m_i.append(round(a, 3))
    return m_i


m_i = full_target(source_data)
#%% Последовательный створ


def sequential_target(source_data):
    m_Δ = ((np.sqrt(2) * source_data['m_B'] * source_data['L'] * 10 ** 3) / 
           (source_data['ro'] * source_data['n+1']))
    
    table_10 = {}
    table_10['(n+1-k)^2'] = []
    table_10['1:(n+1-k)^2'] = []
    table_10['∑1:(n+1-k)^2'] = []
    table_10['sqrt(∑1:(n+1-k)^2)'] = []
    table_10['(n+1-k)*sqrt(∑)'] = []
    table_10['m_δiпр'] = []
    table_10['m_δiоб'] = []
    table_10['m_δi'] = []
    
    flag = 0
    for k in range(source_data['n+1']):
        if k != 0:
            table_10['(n+1-k)^2'].append((source_data['n+1'] - k) ** 2)
            table_10['1:(n+1-k)^2'].append(1 / (source_data['n+1'] - k) ** 2)
            flag += 1 / (source_data['n+1'] - k) ** 2
            table_10['∑1:(n+1-k)^2'].append(flag)
            table_10['sqrt(∑1:(n+1-k)^2)'].append(np.sqrt(flag))
            table_10['(n+1-k)*sqrt(∑)'].append((source_data['n+1'] - k) * np.sqrt(flag))
    
    for k in table_10['(n+1-k)*sqrt(∑)']:
        table_10['m_δiпр'].append(k * m_Δ)
    j = source_data['n'] - 1
    while j > -1:
        table_10['m_δiоб'].append(table_10['m_δiпр'][j])
        j -= 1
    
    for k in range(source_data['n']):
        table_10['m_δi'].append((table_10['m_δiпр'][k] * table_10['m_δiоб'][k]) /
                np.sqrt(table_10['m_δiпр'][k] ** 2 + table_10['m_δiоб'][k] ** 2))
    return table_10, m_Δ


t_10, m_Δ = sequential_target(source_data)
#%% Перекрывающийся створ


def overlapping_target(source_data):
    m_Δ = ((np.sqrt(2) * source_data['m_B'] * source_data['L'] * 10 ** 3) / 
           (source_data['ro'] * source_data['n+1']))
    
    table_11 = {}
    table_11['i:3(n+1)'] = []
    table_11['(i-1)[2i^2-2i(2n+1)-(n+1)]'] = []
    table_11['in(2n+1)'] = []
    table_11['sqrt(A)'] = []
    table_11['m_δ'] = []
    
    for k in range(source_data['n+1']):
        if k != 0:
            table_11['i:3(n+1)'].append(k / (3 * source_data['n+1']))
            table_11['(i-1)[2i^2-2i(2n+1)-(n+1)]'].append((k - 1) * (2 * k ** 2 -
                            2 * k *(2 * source_data['n'] + 1) - source_data['n+1']))
            table_11['in(2n+1)'].append(k * source_data['n'] * (2 * 
                     source_data['n'] + 1))
    for j in range(source_data['n']):
        table_11['sqrt(A)'].append(np.sqrt(table_11['i:3(n+1)'][j] * (
                table_11['(i-1)[2i^2-2i(2n+1)-(n+1)]'][j] + 
                table_11['in(2n+1)'][j])))
        table_11['m_δ'].append(table_11['sqrt(A)'][j] * m_Δ)
        
    return table_11


t_11 = overlapping_target(source_data)
#%% График


mpl.rcParams['font.family'] = 'fantasy'
mpl.rcParams['font.fantasy'] = 'Times New Roman'
x = np.arange(1, source_data['n+1'], 1)
plt.figure(figsize=(15, 7.5)) 
plt.plot(x, m_i, label=r'Схема полного створа')
plt.plot(x, t_10['m_δi'], label=r'Схема последовательных створов')
plt.plot(x, t_11['m_δ'], label=r'Схема перекрывающихся створов')
plt.xlabel(r'$n$', fontsize=20)
plt.ylabel(r'$\mu_{\delta i}$', fontsize=20)
plt.title("""Графики СКО измеренных нестворностей (определённых ошибками визирования)
          для трёх схем построения створа""", fontsize=16)
plt.grid(True)
plt.legend(loc='best', fontsize=13)
plt.savefig('figure_with_legend.png')
plt.show()





























