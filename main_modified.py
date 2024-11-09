#!/usr/bin/env python
# coding: utf-8

import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import matplotlib.pyplot as plt

π = np.pi
def normalize(v):
    return (v / np.linalg.norm(v))

# Matrices de transformations et rotations
def Trans_z(d):
    return np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0,  1, d], [0, 0, 0, 1]], dtype=np.float64)

def Rot_z(θ):
    return np.array([[np.cos(θ), -np.sin(θ), 0, 0], [np.sin(θ), np.cos(θ), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.float64)

def Trans_x(a):
    return np.array([[1, 0, 0, a], [0, 1, 0, 0], [0, 0,  1, 0], [0, 0, 0, 1]], dtype=np.float64)

def Rot_x(α):
    return np.array([[1, 0, 0, 0], [0, np.cos(α), -np.sin(α), 0], [0, np.sin(α), np.cos(α), 0], [0, 0, 0, 1]], dtype=np.float64)

def DH(d, θ, a, α):
    return np.array(Trans_z(d) @ Rot_z(θ) @ Trans_x(a) @ Rot_x(α))

class Chain:
    def __init__(self):
        self.chain = []
        self.base = [0,0,0,1]  # Position de la base
        self.pos = [[0,0]]
        self.target = [1,1]  # Position de la cible
        self.w = 0
        
    # Matrice de transformation de Denavit-Hartenberg
    def chain_DH(self, c):
        return DH(0, c[0], c[1], 0)
    
    # Ajout d'une articulation à notre chaine
    def addArticulatedChain(self, θ, size=1):
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
    
    
    #https://medium.com/unity3danimation/overview-of-jacobian-ik-a33939639ab2
    def getRotationAxis(self, i):
        return [0,0,1]
    
    def getJacobianTranspose(self):
        J = []
        print(len(self.chain))
        for i in range(len(self.chain)):
            # Calcul des coordonnées de l'articulation i
            pos_joint = self.pos[i]
            # Différence de position entre l'articulation i et l'effecteur final
            print("pos_joint: ", pos_joint)
            print("End effector: ", self.getEndEffector())
            delta_pos = np.array(self.getEndEffector()[:2]) - np.array(pos_joint[:2]) #on ne prend que les deux premières valeurs de pos_joint
            # Calcul de la dérivée par rapport à θ_i (2D)
            J_i = np.array([-delta_pos[1], delta_pos[0]])  # Jacobienne 2D pour rotation autour de Z
            J.append(J_i)
        Jt = np.array(J)
        return Jt
    

    def getDeltaOrientation(self):
        end_effector_pos = self.getEndEffector()[:2]  # Prendre seulement x et y
        V = np.array(self.target) - np.array(end_effector_pos)  # Différence entre la cible et l'effecteur
        Jt = self.getJacobianTranspose()  # Jacobienne transposée
        print("V: ", V)
        print("Jt: ", Jt)
        dO = np.dot(Jt, V)  # Calcul de dtheta : variation de l'angle de chaque articulation
        return dO


    #forward inchangée
    def forwardKinematics(self):
        # Calcul de la position de chaque joint
        for i in range(1, len(self.chain)+1): #on itére sur les articulations
            self.pos[i] = self.chain_DH(self.chain[0]) #on initialise la position de l'articulation avec la matrice de transformation de la première articulation
            for j in range(i-1): #on itére sur les articulations précédentes
                # rotation en z d'angle θ, translation en x de size
                self.pos[i] = self.pos[i] @ self.chain_DH(self.chain[j+1]) #on multiplie la matrice de transformation de l'articulation courante par la matrice de transformation de l'articulation précédente
            self.pos[i] = self.pos[i] @ self.base #on multiplie la matrice de transformation de l'articulation courante par la position de la base

    def inverseKinematicsCCD(self, h=0.1):
        dO = self.getDeltaOrientation()
        print("dO: ", dO)
        for i in range(len(self.chain)):
            self.chain[i][0] += dO[i] * h  # Mise à jour de l'angle de chaque articulation
        self.forwardKinematics()  # Recalcule les positions après mise à jour des articulations   

    def plot(self):
        plt.ioff()
        fig = plt.figure(figsize=(10,5))
        plt.xlim(-self.w,self.w)
        plt.ylim(0,self.w)
        plt.xticks([])
        plt.yticks([])
        
        ax = plt.gca()
        for i in range(len(self.pos)-1):
            plt.plot([self.pos[i][0], self.pos[i+1][0]], [self.pos[i][1], self.pos[i+1][1]], zorder=1, linewidth=5.0)
            plt.scatter(self.pos[i][0], self.pos[i][1], color="purple", zorder=2, s=75.0)
        
        # Un cercle pour représenter la cible
        t = plt.Circle((self.target[0], self.target[1]), 0.1, color='r', fill=False)
        ax.add_artist(t)
        return fig


# Exemple d'affichage

chain = Chain()
chain.setStart([0,0])
chain.addArticulatedChain(π/2, 1)
chain.addArticulatedChain(-0.1, 0.5)
chain.addArticulatedChain(-0.1, 1)
chain.addArticulatedChain(-0.1, 0.75)
chain.addArticulatedChain(-0.1, 0.25)
chain.setTarget([-1, 1.0])

ϵ = 0.1 # Erreur souhaitée
max_iter = 100  # Nombre d'itérations maximum
er = 1  # Erreur actuelle
j = 0  # Compteur d'itérations

while (er > ϵ and j < max_iter):
    j = j + 1
    chain.inverseKinematicsCCD()
    chain.forwardKinematics()
    chain.plot()
    f = './res/test' + str(j) + '.png'
    plt.savefig(f)
    ef = chain.getEndEffector()
    dx = chain.target[0] - ef[0]
    dy = chain.target[1] - ef[1]
    er = np.sqrt(dx*dx + dy*dy)
