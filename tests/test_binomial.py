import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import numpy as np 
from scipy.special import binom

from binomial_option_pricer import binomial_option_pricer, call_payoff, put_payoff


# Random function
rng = np.random.default_rng()

tests = 1000
for i in range(tests):

    # Generate random parameters
    S = 100.0+400.0*rng.random()
    K = S+15.0*rng.standard_normal()
    T = rng.integers(1,11)
    N = rng.integers(1,11)
    sigma = rng.random()
    r = rng.random()

    # Calculate these to check no-arbitrage
    R = np.exp(r*T/N)
    u = np.exp(sigma*np.sqrt(T/N))
    d = 1/u
    q = (R-d)/(u-d)

    try:
        A = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r)
        if R <= d or R >= u:
            print("Was expecting no-arbitrage error")
            quit()
    except:
        if R <= d or R >= u:
            continue
        elif R > d and R < u:
            print("Expecting this to pass but didn't")
            quit()

    
    # Multi-period Euro Call test
    C0 = 1.0/R**N * np.sum([(binom(N,j) * q**j * (1-q)**(N-j) * call_payoff(K, S * u**j * d**(N-j))) for j in range(N+1)])
    assert(np.isclose(C0,A.V0))


    # Multi-period American vs Euro Call test
    B = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,european=False)
    assert(A.V0 <= B.V0)

    # C - P = S - K*exp(-rT)
    A_put = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=False)
    assert(np.isclose(A.V0 - A_put.V0, S-K*np.exp(-r*T)))

print(f"All {tests} tests passed!")