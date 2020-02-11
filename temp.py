from __future__ import division
import numpy.matlib
import numpy as np
import math


n_base = 4
u = 1  #taxa de substituição de base por tempo
# priori para a raiz
priori = np.matrix(np.full((n_base, 1), 1.0/n_base))



class Especie:
    def __init__(self, meu_pai, meu_valor, minha_priori = None, meu_tempo = 0):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.num_codon = self.valor.size
        self.L_condicional = np.full((self.num_codon,n_base), None)
        self.P_trs = np.full((self.num_codon,n_base), None)
        if(meu_pai):
          meu_pai.filhos.append(self)
          self.priori = meu_pai.priori
          self.tempo = meu_tempo
          self.trs = self.cria_transicao()
        else:
          self.priori = minha_priori
          
    def cria_transicao(self):
        P = np.identity(n_base)*(math.exp(-u*self.tempo))
        for i in range(P[0,].size):
            P[i,] += ((1 - math.exp(-u*self.tempo))*self.priori[i,0])
        return P.transpose()
    
    def cria_L_condicional_vetor(self):
      self.L_arvore = np.zeros(self.num_codon)
      for i in range(self.num_codon):
        self.cria_L_condicional(i)
    
    def cria_L_condicional(self,i):
      if(not(self.filhos)):# self eh uma folha
          valor = np.zeros(n_base)
          valor[self.valor[i]] = 1
          self.L_condicional[i] = valor
      elif(not(self.pai)): # self eh a raiz
          for filho in self.filhos:
              filho.cria_L_condicional(i)
          self.L_condicional[i] = np.multiply(self.filhos[0].P_trs[i], self.filhos[1].P_trs[i])
          self.P_trs[i] = np.multiply(self.priori.transpose(), self.L_condicional[i])    ##comentar erro com prof
          self.L_arvore[i] = np.sum(self.P_trs[i]) 
      else: # self eh no interno
          for filho in self.filhos:
              filho.cria_L_condicional(i)
              filho.P_trs[i] = np.dot(filho.L_condicional[i],filho.trs).transpose()
          self.L_condicional[i] = np.multiply(self.filhos[0].P_trs[i], self.filhos[1].P_trs[i])
          self.P_trs[i] = np.dot(self.L_condicional[i], self.trs).transpose()



#S0 = Especie(None, np.full(2,None), minha_priori = priori)
#6 = Especie(S0, np.full(2,None), meu_tempo = 6)
#S1 = Especie(S6, np.array([1,2]), meu_tempo = 1)
#S2 = Especie(S6, np.array([0,2]), meu_tempo = 0.7 )
#S8 = Especie(S0, np.full(2,None), meu_tempo = 0.4)
#S3 = Especie(S8, np.array([0,3]), meu_tempo = 0.6)
#S7 = Especie(S8, np.full(2,None), meu_tempo = 0.9)
#S4 = Especie(S7, np.array([0,2]), meu_tempo = 0.4)
#S5 = Especie(S7, np.array([1,3]), meu_tempo = 0.3)
#S0.cria_L_condicional_vetor()
#print(S1.P_trs)
##print(S6.L_condicional)
#print(S0.L_condicional)
#print(S0.L_arvore)

h = np.array([True, False, False, False])
A = np.where(h)[0]
print(A)
b  = np.array([[1,2],[3,4],[5,6],[6,5],[1,0]])

a = np.array([[0.1137,0.0342,0.0555,0.0342]])
h = np.array([[0.4173,0.1942,0.1942,0.1942],
[0.1942,0.4173,0.1942,0.1942],
[0.1942,0.1942,0.4173,0.1942],
[0.1942,0.1942,0.1942,0.4173]])
print(np.dot(a,h))