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
    assert(round(A.V0,5) <= round(B.V0,5))

# Intrinsic value tests; price of American option should be at least the payoff of immediate exercising of option
def test_intrinsic_value(S,K,T,N,sigma,r):
    american_call = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=True,european=False)
    american_put = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=False,european=False)
    european_call = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=True,european=True)
    european_put = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=False,european=True)

    assert(american_call.V0 > max(S-K,0) or np.isclose(american_call.V0,max(S-K,0)))
    assert(american_put.V0 > max(K-S,0) or np.isclose(american_put.V0,max(K-S,0)))
    assert(european_call.V0 > max(S-K*np.exp(-r*T),0) or np.isclose(european_call.V0,max(S-K*np.exp(-r*T),0)))
    assert(european_put.V0 > max(K*np.exp(-r*T)-S,0) or np.isclose(european_put.V0,max(K*np.exp(-r*T)-S,0)))

# Test no-arbitrage price bounds for European options
def test_no_arbitrage_price_bounds(S,K,T,N,sigma,r):
    C = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=True)
    P = binomial_option_pricer(S=S,K=K,T=T,N=N,sigma=sigma,r=r,call=False)

    assert(C.V0 < S or np.isclose(C.V0,S))
    assert(P.V0 < K*np.exp(-r*T) or np.isclose(P.V0, K*np.exp(-r*T)))

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

    test_intrinsic_value(S,K,T,N,sigma,r)

    if european: test_no_arbitrage_price_bounds(S,K,T,N,sigma,r)

print(f"All {TESTS} tests passed!")