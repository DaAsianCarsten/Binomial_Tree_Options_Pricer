import numpy as np 
from scipy.stats import norm

def black_scholes_european_price(S,K,r,T,sigma,call=True):
    d1 = (np.log(S/K) + (r+sigma**2/2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)

    # For a call
    V = S*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)

    # For a put; use C-P = S-K*exp(-rT)
    if not call:
        V = V - (S - K*np.exp(-r*T))

    return V