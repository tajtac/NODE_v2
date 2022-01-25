import numpy as np

class neoHook(): #Incompressible neo-Hookean
    def __init__(self, param):
        self.C1 = param[0]
    def kinematics(self, lm):
        n = np.size(lm,0)
        F = np.zeros([n,3,3])
        F[:,0,0] = lm[:,0]
        F[:,1,1] = lm[:,1]
        F[:,2,2] = 1/(lm[:,0]*lm[:,1])
        C = F*F #Since F is diagonal
        I1 = C.trace(axis1=1, axis2=2)
        return F, C, I1
    def Psi(self, lm):
        C1 = self.C1
        _,_,I1 = self.kinematics(lm)
        return C1*(I1-3)
    def S(self, lm):
        C1 = self.C1
        _, C, _ = self.kinematics(lm)
        C_inv = np.linalg.inv(C)
        n = C.shape[0]
        I = np.array([np.identity(3) for i in range(n)])
        p = 2*C1*C[:,2,2]
        S = 2*C1*I - p[:, None, None]*C_inv
        return S
    def sigma(self, lm):
        F, _, _, = self.kinematics(lm)
        S = self.S(lm)
        sigma = np.einsum('...ij,...jk,...lk->...il', F, S, F)
        return sigma
        

class MR(): #Source: Continuummechanics
#Assumptions: Fully incompressible material. Thus paramter D is irrelevant.
    #Strain Energy
    def __init__(self, params):
        self.C10 = params[0]
        self.C01 = params[1]
        self.C20 = params[2]
    def Psi(self, lm): #lm.shape = (n,2)
        C10 = self.C10
        C01 = self.C01
        C20 = self.C20
        lm3 = np.zeros(lm.shape[0])
        lm3[:] = 1/(lm[:,0] * lm[:,1])
        I1 = lm[:,0]**2 + lm[:,1]**2 + lm3**2
        I2 = 1/lm[:,0]**2 + 1/lm[:,1]**2 + 1/lm3**2
        return C10*(I1-3) + C01*(I2-3) + C20*(I1-3)**2
    #stress tensor given lm1 and lm2, assuming sigma3=0 and J=1
    def sigma(self, lm): #lm.shape = (n,2)
        C10 = self.C10
        C01 = self.C01
        C20 = self.C20
        lm1 = lm[:,0]
        lm2 = lm[:,1]
        lm3 = 1/(lm1*lm2)
        sigma1 = 2*(C10*(lm1**2 - lm3**2) - C01*(1/lm1**2 - 1/lm3**2) +
                  2*C20*(lm1**2 - lm3**2)*(lm1**2 + lm2**2 + lm3**2 - 3))
        sigma2 = 2*(C10*(lm2**2 - lm3**2) - C01*(1/lm2**2 - 1/lm3**2) +
                  2*C20*(lm2**2 - lm3**2)*(lm1**2 + lm2**2 + lm3**2 - 3))
        sigma = np.zeros((lm.shape[0],3,3))
        sigma[:,0,0] = sigma1
        sigma[:,1,1] = sigma2
        return sigma

class GOH():
    #Paper: Propagation of material behavior uncertainty in a nonlinear finite
    #element model of reconstructive surgery

    #This assumes fully incompressible material. Therefore bulk modulus, K and 
    #U_vol are 0.
    def __init__(self, params): #lm.shape = (n,2)
        self.mu = params[0]
        self.k1 = params[1]
        self.k2 = params[2]
        self.kappa = params[3]
        self.theta = params[4]
    def kinematics(self, lm):
        theta = self.theta
        n = np.size(lm,0)
        F = np.zeros([n,3,3])
        F[:,0,0] = lm[:,0]
        F[:,1,1] = lm[:,1]
        F[:,2,2] = 1/(lm[:,0]*lm[:,1])
        C = F*F
        I1 = C.trace(axis1=1, axis2=2)
        e_0 = [np.cos(theta), np.sin(theta), 0]
        I4 = np.einsum('i,pij,j->p', e_0, C, e_0)
        return F, C, I1, I4, e_0
    #Strain Energy
    def Psi(self, lm): #lm.shape = (n,2)
        mu = self.mu
        k1 = self.k1
        k2 = self.k2
        kappa = self.kappa
        _, _, I1, I4, _ = self.kinematics(lm)
        E = kappa*(I1-3) + (1-3*kappa)*(I4-1)
        Psi_iso = mu/2*(I1-3)
        #U_vol = 0.0 because we are assuming fully incompressible material
        Psi_aniso = k1/2/k2*(np.exp(k2*E**2) - 1)
        Psi = Psi_iso + Psi_aniso
        return Psi
    #Partial derivatives of Strain Energy wrt Invariants I1, I4
    def partial(self, lm):
        mu = self.mu
        k1 = self.k1
        k2 = self.k2
        kappa = self.kappa
        _, _, I1, I4, _ = self.kinematics(lm)
        E = kappa*(I1-3) + (1-3*kappa)*(I4-1)

        Psi1 = mu/2 + k1*np.exp(k2*E**2)*E*kappa
        Psi4 = k1*np.exp(k2*E**2)*E*(1-3*kappa)
        return Psi1, Psi4
    def S(self, lm): 
        mu = self.mu
        k1 = self.k1
        k2 = self.k2
        kappa = self.kappa
        F, C, I1, I4, e_0 = self.kinematics(lm)
        eiej = np.outer(e_0,e_0) #e_i dyadic e_j
        E = kappa*(I1-3) + (1-3*kappa)*(I4-1)
        C_inv = np.linalg.inv(C)
        n = np.size(lm,0)
        I = np.identity(3)
        S_iso = np.zeros([n,3,3])
        S_vol = np.zeros([n,3,3])
        for i in range(0,3):
            for j in range(0,3):
                S_iso[:,i,j] = mu*(I[i,j] - 1/3*I1[:]*C_inv[:,i,j])
        #There is a mistake in the paper in the equation above.
        #The J**(-2/3) is supposed to be factored out. See Prof. Tepole's notes.
        #S_vol = 0.0 because we are assuming fully incompressible material
        eiej = np.outer(e_0,e_0) #e_i dyadic e_j
        dI1dC = np.zeros(n)
        dI4dC = np.zeros(n)
        S_aniso = np.zeros([n,3,3])
        for i in range(0,3):
            for j in range(0,3):
                dI1dC[:] = I[i,j] - 1/3*I1[:]*C_inv[:,i,j]
                # dI1dC[:] = I[i,j] #Both this and the line above are true when C_isoch = C
                dI4dC[:] = eiej[i,j] - 1/3*I4[:]*C_inv[:,i,j]
                # dI4dC[:] = eiej[i,j] #Both this and the line above are true when C_isoch = C

                S_aniso[:,i,j] = 2*k1*np.exp(k2*E[:]**2)*E[:]*(kappa*dI1dC[:] + (1-3*kappa)*dI4dC[:])
        p = -(S_iso[:,2,2] + S_aniso[:,2,2])*C[:,2,2]
        for i in range(0,3):
            for j in range(0,3):
                S_vol[:,i,j] = p[:]*C_inv[:,i,j]
        S = S_iso + S_aniso + S_vol
        return S
    def sigma(self, lm):
        F, _, _, _, _ = self.kinematics(lm)
        S = self.S(lm)
        sigma = np.einsum('...ij,...jk,...lk->...il', F, S, F)
        return sigma
    def S_from_F(self, F): #F.shape=[n,3,3]
        mu = self.mu
        k1 = self.k1
        k2 = self.k2
        kappa = self.kappa
        theta = self.theta
        C = np.einsum('...ji,...jk->...ik', F, F)
        I1 = C.trace(axis1=1, axis2=2)
        e_0 = [np.cos(theta), np.sin(theta), 0]
        I4 = np.einsum('i,pij,j->p', e_0, C, e_0)
        eiej = np.outer(e_0,e_0) #e_i dyadic e_j
        E = kappa*(I1-3) + (1-3*kappa)*(I4-1)
        C_inv = np.linalg.inv(C)
        n = F.shape[0]
        I = np.identity(3)
        S_iso = np.zeros([n,3,3])
        S_vol = np.zeros([n,3,3])
        for i in range(0,3):
            for j in range(0,3):
                S_iso[:,i,j] = mu*(I[i,j] - 1/3*I1[:]*C_inv[:,i,j])
        #There is a mistake in the paper in the equation above.
        #The J**(-2/3) is supposed to be factored out. See Prof. Tepole's notes.
        #S_vol = 0.0 because we are assuming fully incompressible material
        eiej = np.outer(e_0,e_0) #e_i dyadic e_j
        dI1dC = np.zeros(n)
        dI4dC = np.zeros(n)
        S_aniso = np.zeros([n,3,3])
        for i in range(0,3):
            for j in range(0,3):
                dI1dC[:] = I[i,j] - 1/3*I1[:]*C_inv[:,i,j]
                # dI1dC[:] = I[i,j] #Both this and the line above are true when C_isoch = C
                dI4dC[:] = eiej[i,j] - 1/3*I4[:]*C_inv[:,i,j]
                # dI4dC[:] = eiej[i,j] #Both this and the line above are true when C_isoch = C

                S_aniso[:,i,j] = 2*k1*np.exp(k2*E[:]**2)*E[:]*(kappa*dI1dC[:] + (1-3*kappa)*dI4dC[:])
        p = -(S_iso[:,2,2] + S_aniso[:,2,2])/C_inv[:,2,2]
        for i in range(0,3):
            for j in range(0,3):
                S_vol[:,i,j] = p[:]*C_inv[:,i,j]
        S = S_iso + S_aniso + S_vol
        return S

class HGO():
    #Ref: M. Liu et al 2020
    def __init__(self, params):
        self.C10   = params[0]
        self.k1    = params[1]
        self.k2    = params[2]
        self.theta = params[3]
    def kinematics(self, lm):
        theta = self.theta
        v0 = np.array([np.cos(theta), np.sin(theta), 0])
        w0 = np.array([np.cos(theta),-np.sin(theta), 0])
        V0 = np.outer(v0,v0)
        W0 = np.outer(w0,w0)

        n = lm.shape[0]
        F = np.zeros([n,3,3])
        F[:,0,0] = lm[:,0]
        F[:,1,1] = lm[:,1]
        F[:,2,2] = 1/(lm[:,0]*lm[:,1])
        C = F*F #Since F^T = F.
        I1 = np.trace(C)
        I4 = np.tensordot(C,V0)
        I6 = np.tensordot(C,W0)

        return F, C, I1, I4, I6, V0, W0
    def S(self, lm):
        C10   = self.C10
        k1    = self.k1
        k2    = self.k2
        theta = self.theta
        F, C, I1, I4, I6, V0, W0 = self.kinematics(lm)

        invC = np.linalg.inv(C)
        I = np.eye(3)
        Psi1 = C10
        Psi4 = k1*(I4-1)*np.exp(k2*(I4-1)**2)
        Psi6 = k1*(I6-1)*np.exp(k2*(I6-1)**2)
        
        S2 = 2*Psi1*I + 2*Psi4[:, None, None]*V0 + 2*Psi6[:, None, None]*W0
        p = S2[:,2,2]/invC[:,2,2]
        S = S2 - p[:, None, None]*invC
        return S
    def sigma(self, lm):
        F, _, _, _, _, _, _ = self.kinematics(lm)
        S = self.S(lm)
        sigma = np.einsum('...ij,...jk,...lk->...il', F, S, F)
        return sigma
    def Psi(self, lm):
        C10   = self.C10
        k1    = self.k1
        k2    = self.k2
        theta = self.theta
        _, _, I1, I4, I6, _, _ = self.kinematics(lm)

        invC = np.linalg.inv(C)
        I = np.eye(3)

        Psi = C10*(I1-3) + k1/2/k2*(np.exp(k2*(I4-1)**2)-1 + np.exp(k2*(I6-1)**2)-1)
        return Psi
        
    
class Fung():
    #Source: Fung et al. 1979. Replace F with Q and C with c1 to avoid confusion with deformation tensors.
    #Also, set * values to zero.
    def __init__(self, params):
        self.c1 = params[0]
        self.a1 = params[1]
        self.a2 = params[2]
        self.a4 = params[3]
    def kinematics(self, lm):
        n = lm.shape[0]
        F = np.zeros([n,3,3])
        F[:,0,0] = lm[:,0]
        F[:,1,1] = lm[:,1]
        F[:,2,2] = 1/(lm[:,0]*lm[:,1])
        C = F*F
        E_11 = 0.5*(C[:,0,0]-1)
        E_22 = 0.5*(C[:,1,1]-1)
        E_Z = E_11
        E_theta = E_22
        return F, E_Z, E_theta

    def S(self, lm):
        c1 = self.c1
        a1 = self.a1
        a2 = self.a2
        a4 = self.a4
        F, E_Z, E_theta = self.kinematics(lm)

        Q = a1*E_theta**2 + a2*E_Z**2 + 2*a4*E_theta*E_Z
        S_theta = c1*(a1*E_theta + a4*E_Z)*np.exp(Q) #Eq. (4)
        S_Z     = c1*(a4*E_theta + a2*E_Z)*np.exp(Q) #Eq. (4)
        S = np.zeros((lm.shape[0],3,3))
        S[:,0,0] = S_theta
        S[:,1,1] = S_Z
        return S
    
    def sigma(self, lm):
        F, _, _ = self.kinematics(lm)
        S = self.S(lm)
        sigma = np.einsum('...ij,...jk,...lk->...il', F, S, F)
        return sigma
        
    def Psi(self, lm):
        c1 = self.c1
        a1 = self.a1
        a2 = self.a2
        a4 = self.a4
        _, E_Z, E_theta = self.kinematics(lm)

        Q = a1*E_theta**2 + a2*E_Z**2 + 2*a4*E_theta*E_Z
        Psi = c1/2*(np.exp(Q)-1)
        return Psi




