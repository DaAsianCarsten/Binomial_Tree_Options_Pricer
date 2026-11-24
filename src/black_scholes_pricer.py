import numpy as np 
from scipy.stats import norm

def black_scholes_price(S,K,r,T,sigma,q_d=0,call=True,european=True):
    if (q_d != 0 and not european): raise TypeError("This model doe not apply to American-style options with nonzero dividend yield due to early exercising.")
    if (not european and not call): raise TypeError("This model does not apply to American-style put options due to early exercising.")

    d1 = (np.log(S/K) + (r-q_d+sigma**2/2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)

    # For a call
    if call: V = S*np.exp(-q_d*T)*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
    # For a put
    else: V = K*np.exp(-r*T)*norm.cdf(-d2) - S*np.exp(-q_d*T)*norm.cdf(-d1)

    return V