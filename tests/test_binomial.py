import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import numpy as np 
from scipy.special import binom

from binomial_option_pricer import binomial_option_pricer, call_payoff, put_payoff
from black_scholes_pricer import black_scholes_european_price

# Random function
rng = np.random.default_rng()

# Number of performed tests
TESTS = 1000

# Testing the fact for European options, C-P = S - K*exp(-rT)
def test_put_call_parity(S,K,T,N,sigma,r,call):
    A = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=True)
    B = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=False)

    assert(np.isclose(A.V0 - B.V0, S-K*np.exp(-r*T)))

# Test that the no-arbitrage rule flag gets raised when appropriate
def test_no_arbitrage(S,K,T,N,sigma,r,call):
    R = np.exp(r*T/N)
    u = np.exp(sigma*np.sqrt(T/N))
    d = 1/u
    q = (R-d)/(u-d)

    rule_followed = (0<q<1)
    test_rule_followed = True
    try:
        A = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r)
    except:
        test_rule_followed = False

    assert(rule_followed == test_rule_followed)
    return rule_followed

# Compare the European price in model follows the analytical discrete-time calculation
def multi_period_euro_test(S,K,T,N,sigma,r,call):
    R = np.exp(r*T/N)
    u = np.exp(sigma*np.sqrt(T/N))
    d = 1/u
    q = (R-d)/(u-d)

    if call: g = call_payoff
    else: g = put_payoff

    A = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=call)
    C0 = 1.0/R**N * np.sum([(binom(N,j) * q**j * (1-q)**(N-j) * g(K, S * u**j * d**(N-j))) for j in range(N+1)])
    assert(np.isclose(C0,A.V0))

# Compare American vs European style options with exact parameters; American option should NOT be priced lower
def test_american_european_comparison(S,K,T,N,sigma,r,call):
    A = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=call,european=True)
    B = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=call, european=False)
    assert(A.V0 <= B.V0)


for i in range(TESTS):

    # Generate random parameters
    S = 1.0+1000.0*rng.random()
    K = 1.0+1000.0*rng.random()
    T = rng.integers(1,101)
    N = rng.integers(1,101)
    sigma = rng.random()
    r = rng.random()
    call = rng.integers(0,2)
    european = rng.integers(0,2)

    # Perform appropriate tests
    b = test_no_arbitrage(S,K,T,N,sigma,r,call)
    if not b: continue

    if european: test_put_call_parity(S,K,T,N,sigma,r,call)

    if european: multi_period_euro_test(S,K,T,N,sigma,r,call)

    test_american_european_comparison(S,K,T,N,sigma,r,call)

print(f"All {TESTS} tests passed!")