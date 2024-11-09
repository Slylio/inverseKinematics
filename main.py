#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
np.seterr(divide='ignore', invalid='ignore')
import matplotlib.pyplot as plt


# In[2]:


π = np.pi
def normalize(v):
    return (v / np.linalg.norm(v))


# ![title](images/TZDZ.png)

# In[3]:

# Matrices de transformations et rotations
def Trans_z(d):
    return np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0,  1, d], [0, 0, 0, 1]], dtype=np.float64)

def Rot_z(θ):
    return np.array([[np.cos(θ), -np.sin(θ), 0, 0], [np.sin(θ), np.cos(θ), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float64)


# ![title](images/TXDX.png)

# In[4]:

# Matrices de transformations et rotations
def Trans_x(a):
    return np.array([[1, 0, 0, a], [0, 1, 0, 0], [0, 0,  1, 0], [0, 0, 0, 1]], dtype=np.float64)

def Rot_x(α):
    return np.array([[1, 0, 0, 0], [0, np.cos(α), -np.sin(α), 0], [0, np.sin(α),  np.cos(α), 0], [0, 0, 0, 1]], dtype=np.float64)


# ![title](images/DH.png)

# In[5]:

"""
Fonction de transformation de Denavit-Hartenberg : permet de calculer la matrice de transformation
d : distance sur l'axe z
θ : angle de rotation sur l'axe z
a : distance sur l'axe x
α : angle de rotation sur l'axe x
"""
def DH(d, θ, a, α):
    return np.array(Trans_z(d) @ Rot_z(θ) @ Trans_x(a) @ Rot_x(α))


# # CCD
# Pour chaque joint du end effector à la base :<br>
# Avec p la position du joint courrant, e la position du end effector et t la position de la cible.<br>
# <div style="align:left">Soit α l'angle entre $\vec{pe}$ et $\vec{pt}$</div>
# On cherche l'angle α tel que :
# $
# \cos \alpha = \frac{ (\mathbf{e} - \mathbf{p}) \cdot (\mathbf{t} - \mathbf{p})}{|\mathbf{e} - \mathbf{p}| |\mathbf{t} - \mathbf{j}| }
# $
# et
# $
# sin \alpha = \frac {(e_x - p_x)(t_y - p_y)-(e_y - p_y)(t_x - p_x)} {|\mathbf{e} - \mathbf{p}| |\mathbf{t} - \mathbf{p}|}
# $

# In[8]:

# Classe pour la chaine articulée
class Chain:
    def __init__(self):
        self.chain = []
        self.base = [0,0,0,1] # Position de la base
        self.pos = [[0,0]]
        self.target = [1,1] # Position de la cible
        self.w = 0
        
    # Matrice de transformation de Denavit-Hartenberg avec d=0 et α=0 car on rota sur z et on translat sur x
    def chain_DH(self, c):
        return DH(0, c[0], c[1], 0)
    
    """
    Ajout d'une articulation à notre chaine
    θ : angle de rotation
    size : taille de l'articulation
    """
    def addArticulatedChain(self, θ, size=1) :
        self.chain.append([θ, size])
        self.pos.append([0,0]) 
        self.w += size
    
    # Position de départ
    def setStart(self, p):
        self.pos[0] = p
        
    # Position de l'end effector
    def getEndEffector(self):
        return self.pos[-1]
    
    # Position de la cible
    def setTarget(self, target):
        self.target = target
        
    """
    Calcul de la cinématique direct

    Pour chaque articulation, on calcule la matrice de transformation de Denavit-Hartenberg, on la multiplie par la matrice de transformation de l'articulation précédente et on l'applique à la position de la base
    """
    def forwardKinematics(self):
        for i in range(1, len(self.chain)+1): #on itére sur les articulations
            self.pos[i] = self.chain_DH(self.chain[0]) #on initialise la position de l'articulation avec la matrice de transformation de la première articulation
            for j in range(i-1): #on itére sur les articulations précédentes
                # rotation en z d'angle θ, translation en x de size
                self.pos[i] = self.pos[i] @ self.chain_DH(self.chain[j+1]) #on multiplie la matrice de transformation de l'articulation courante par la matrice de transformation de l'articulation précédente
            self.pos[i] = self.pos[i] @ self.base #on multiplie la matrice de transformation de l'articulation courante par la position de la base
            
    """
    Calcul de la cinématique inverse
    
    clamp : limite de l'angle maximal en radian
    """
    def inverseKinematicsCCD(self, clamp=0.1):
        #on prend la position de l'end effector
        e = np.array(self.pos[-1])

        #On itère de l'extrémité à la base
        for i in range(len(self.pos)-1,-1,-1): 
            # Articulation courante
            p = np.array(self.pos[i])
            
            # Calcul des vecteurs entre l'articulation courante et l'end effector et entre l'articulation et la cible
            c = e[0:2] - p[0:2]
            t = self.target - p[0:2]
            
            #On normalise les vecteurs
            #On ajoute un 0 pour la dimension z (pouvoir faire des calculs 3D dans un plan 2D)
            c = normalize(np.append(c,0))
            t = normalize(np.append(t,0))
            
            #
            cθ = np.dot(c,t)
            if cθ < 1.0 :# 1.0 -> θ = 0
                # Cross for SIN
                k = np.cross(t,c)
                # On limite θ à une valeur maximal δ
                θ = np.minimum(clamp, np.arccos(cθ))
                if k[2] > 0 :
                    # Angle positive = anti-horaire
                    self.chain[i][0] = self.chain[i][0] - θ
                else :
                    # Angle negatif = horaire
                    self.chain[i][0] = self.chain[i][0] + θ
        
    def plot(self) :
        plt.ioff()
        fig = plt.figure(figsize=(10,5))
        plt.xlim(-self.w,self.w)
        # 0 For the tests, ca be change to w/e we want
        plt.ylim(0,self.w)
        # remove to see coords
        plt.xticks([])
        plt.yticks([])
        
        ax = plt.gca()

        for i in range(len(self.pos)-1):
            plt.plot([self.pos[i][0], self.pos[i+1][0]], [self.pos[i][1], self.pos[i+1][1]],zorder=1, linewidth = 5.0)
            plt.scatter(self.pos[i][0], self.pos[i][1], color="purple", zorder=2, s=75.0)
                

        # Un cercle pour representer la cible
        t = plt.Circle((self.target[0], self.target[1]), 0.1, color='r', fill=False)
        ax.add_artist(t)
        return fig


# ## Exemple d'affichage

# In[9]:


chain = Chain()
# Where the arm starts
chain.setStart([0,0])
# First Joint
chain.addArticulatedChain(π/2, 1)
chain.addArticulatedChain(-0.1, 0.5)
# Second Joint
chain.addArticulatedChain(-0.1, 1)
# Third Joint
chain.addArticulatedChain(-0.1, 0.75)
# End effector
chain.addArticulatedChain(-0.1, 0.25)
# position de la cible
chain.setTarget([-1, 1.0])

# Conditions d'arrets
# Erreur souhaitée
ϵ = 0.1
# Nombre d'iterations maximum
max_iter = 100

# Current error
er = 1.0
# For image saving
j = 0

while (er > ϵ and j < max_iter):
    j = j + 1
    chain.inverseKinematicsCCD()
    chain.forwardKinematics()
    chain.plot()
    f = './res/test'+str(j)+'.png'
    plt.savefig(f)
    # Erreur = Distance euclidienne entre l'end effector et la target
    ef = chain.getEndEffector()
    dx = chain.target[0] - ef[0]
    dy = chain.target[1] - ef[1]
    er = np.sqrt(dx*dx + dy*dy)








# %%

# %%
