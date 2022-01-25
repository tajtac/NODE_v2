import jax.numpy as jnp
import jax
from jax import grad, vmap, jit, partial, random, jacrev
from jax.experimental.ode import odeint
from jax.experimental import optimizers
from jax.scipy.optimize import minimize
from jax.lax import scan

@jit
def forward_pass(H, Ws):
    N_layers = len(Ws)
    for i in range(N_layers - 1):
        H = jnp.matmul(H, Ws[i])
        H = jnp.tanh(H)
    Y = jnp.matmul(H, Ws[-1])
    return Y

@jit
def NN_old(y0, params):
    f = lambda y, t: forward_pass(jnp.array([y]),params) # fake time argument for ODEint
    return odeint(f, y0, jnp.array([0.0,1.0]))[-1] # integrate between 0 and 1 and return the results at 1

#The same function as NN_old except using Euler integration
@jit
def NN(y0, params, steps = 10):
    body_func = lambda y, i: (y + forward_pass(jnp.array([y]), params)[0], None)
    out, _ = scan(body_func, y0, None, length = steps)
    return out
NN_vmap = vmap(NN, in_axes=(0, None), out_axes=0)

@jit
def NODE_sigma(F, params):
    C = jnp.dot(F.T, F)
    S = NODE_S(C, params)
    return jnp.einsum('ij,jk,kl->il', F, S, F.T)
NODE_sigma_vmap = vmap(NODE_sigma, in_axes=(0, None), out_axes=0)

@jit
def NODE_S(C, params):
    Psi1_params, Psi2_params, Psi4v_params, Psi4w_params, J1_params, J2_params, J3_params, J4_params, J5_params, J6_params, J_weights, theta = params
    w1, w2, w3, w4, w5, w6 = jnp.abs(J_weights)
    v0 = jnp.array([ jnp.cos(theta), jnp.sin(theta), 0])
    w0 = jnp.array([-jnp.sin(theta), jnp.cos(theta), 0])
    V0 = jnp.outer(v0, v0)
    W0 = jnp.outer(w0, w0)
    I1 = jnp.trace(C)
    C2 = jnp.einsum('ij,jk->ik', C, C)
    I2 = 0.5*(I1**2 - jnp.trace(C2))
    I4v = jnp.einsum('ij,ij',C,V0)
    I4w = jnp.einsum('ij,ij',C,W0)
    
    I4v = jnp.max([I4v,0.0])
    I4w = jnp.max([I4w,0.0])
    Cinv = jnp.linalg.inv(C)

    I1 = I1-3
    I2 = I2-3
    Iv = Iv-1
    Iw = Iw-1
    J1 = I1+I2
    J2 = I1+Iv
    J3 = I1+Iw
    J4 = I2+Iv
    J5 = I2+Iw
    J6 = Iv+Iw

    Psi1 = NN(I1,  W1_params)
    Psi2 = NN(I2,  W2_params)
    Psiv = NN(Iv,  Wv_params)
    Psiw = NN(Iw,  Ww_params)
    Phi1 = NN(J1,  J1_params)
    Phi2 = NN(J2,  J2_params)
    Phi3 = NN(J3,  J3_params)
    Phi4 = NN(J4,  J4_params)
    Phi5 = NN(J5,  J5_params)
    Phi6 = NN(J6,  J6_params)
    
    Psiv = np.max([Psiv, 0])
    Psiw = np.max([Psiw, 0])
    Phi1 = np.max([Phi1, 0])
    Phi2 = np.max([Phi2, 0])
    Phi3 = np.max([Phi3, 0])
    Phi4 = np.max([Phi4, 0])
    Phi5 = np.max([Phi5, 0])
    Phi6 = np.max([Phi6, 0])
    
    Psi1 = Psi1 + w1*Phi1 + w2*Phi2 + w3*Phi3
    Psi2 = Psi2 + w1*Phi1 + w4*Phi4 + w5*Phi5
    Psiv = Psiv + w2*Phi2 + w4*Phi4 + w6*Phi6
    Psiw = Psiw + w3*Phi3 + w5*Phi5 + w6*Phi6
    
    p = -C[2,2]*(2*Psi1 + 2*Psi2*(I1 - C[2,2]) + 2*Psiv*V0[2,2] + 2*Psiw*W0[2,2])
    S = p*Cinv + 2*Psi1*jnp.eye(3) + 2*Psi2*(I1*jnp.eye(3)-C) + 2*Psiv*V0 + 2*Psiw*W0
    return S
NODE_S_vmap = vmap(NODE_S, in_axes=0, out_axes=0)