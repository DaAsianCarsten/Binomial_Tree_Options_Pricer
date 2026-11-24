# Binomial_Tree_Options_Pricer
<!-- WRITE OUT THE OVERVIEW MORE!!!!!!!!!!!!!!!!! A discrete-time model that calculates the fair price of an option using risk-neutral valuation. 

## What an option is
    european vs american
    call vs put
    payoff vs profit
    moneyness

## Basic principles on pricing options and (financial) assumptions
    time-value of money
    we assume risk-free rate is known and constant
    assume frictionless trading
    assume assets can be traded/shorted at continuous values w/o restrictions
    no-arbitrage / consistent pricing.
        we assume that market is efficient enough so arbitrage opportunities don't persist
    assume no dividends for now

## Core model assumptions
    price moves like GBM (talk more about it later)
    talk about the CRR parametrization and why it makes sense under our assumptions
    constant volatility throughout option lifespan (of course false) maybe mention stuff about implied volatility and other stochastic vol models

## The binomial model and tree
    the math behind it/breifly talk about CRR, explain later
        talk about d<R<u
        derive risk neutral prob q
    backwards induction/ mention about how its most general; applies to american options
    the algorithm
    talk briefly about observations and just refer to ipynb (breif summary of findings)
    time/space complexity

## Stochastic calculus
    Go over the GBM SDE and solution
        go briefly about Ito's lemma
    Talk about changes in measure/ Girsanov's Theorem

## BS derivation
    derive BS

## Numerical results from BS convergence experiments
    talk about the oscillations and irregular beating wieh S_0 \neq K

## The Greeks
    do very brief defs for all 5 greeks as most are in the ipynb

## Incorporating dividends
    talk about the change in risk netrual prob q, but all aspects of models stay the same
    briefly talk about option price inc/dec for put/call respectively but refer to ipynb (brief summary of findings)

## Model limitations, future work
    volatility is of course not constant, and it isn't consistent i.e. implied volatility surface
    stock price does not move like GBM
    assumed constant market interest rates r (not true)
    assumed frictionless market (not true)

## References
    refs
-->
## *Why American-style Calls Never Benefit from Early Exercising (assuming no dividends)

From the tree model, suppose that one node does benefit from early exercising. For American-style call options, the value of the option at this node is given by:

$$
V=\max\left(g(S),\frac{1}{R}\left(q\cdot g(Su)+(1-q)\cdot g(Sd)\right)\right)
$$

where $g(S)=\max(S-K,0)$ is the payoff of the call option at stock price $S$. Recall that $q$ (the risk-neutral probability) in this case is defined by $q:=\frac{R-d}{u-d}$. We have several cases for where the prices are positioned relative to the strike price $K$.

### Case 1: $K\leq Sd< S < Su $

In this case, all the payoff functions are in the money, so


$$
\frac{1}{R}\left(q\cdot g(Su)+(1-q)\cdot g(Sd)\right) = \frac{1}{R}\left(q(Su-K)+(1-q)(Sd-K)\right)=\frac{1}{R}(RS-K)=S-\frac{K}{R} \qquad (1)
$$


and since $R \geq 1$, we have that $S-\frac{K}{R}\geq S-K=g(S)$. So, exercising early isn't optimal in this case.

### Case 2: $Sd< K\leq S < Su$

In this case, we have $g(Sd)=0$, so

$$
\frac{1}{R}\left(q\cdot g(Su)+(1-q)\cdot g(Sd)\right)=\frac{1}{R}\cdot q(Su-K)\qquad (2)
$$

Observe that

$$
\begin{align*}
\frac{q}{R} &< d\\
\Leftrightarrow \frac{R-d}{Ru-Rd} &< d && (\text{Definition of } q)\\
\Leftrightarrow R-d &< Rud-Rd^2\\
\Leftrightarrow R-d &< R-Rd^2 && (ud=1)\\
\Leftrightarrow -d &< -Rd^2\\
\Leftrightarrow \frac{1}{d} &> R\\
\Leftrightarrow R &< u && (u=1/d)\\
\end{align*}
$$

which is true because of our no-arbitrage assumption that $d<R<u$. So, plugging in the result that $\frac{q}{R}<d$ into (2), we get that

$$
\frac{1}{R}\cdot q(Su-K) < d(Su-K)=S-Kd \qquad (3)
$$

where that last step uses the fact $ud=1$. Since $d < 1$, we have that $S-Kd>S-K=g(S)$. So, in this case, we have that $\frac{1}{R}\left(q\cdot g(Su)+(1-q)\cdot g(Sd)\right) > g(S)$ meaning that exercising early isn't optimal in this case.

### Case 3: $Sd < S < K\leq Su$

In this case, $g(S)=0$, so the comparison is trivial:

$$
\frac{1}{R}\left(q\cdot g(Su)+(1-q)\cdot g(Sd)\right)=\frac{q}{R}(Su-K)>0=g(S)\qquad (4)
$$

So, in this case, early exercising is not optimal.

### Case 4: $Sd < S < Su < K$

In this case, $g(S)=g(Su)=g(Sd)=0$, making the comparison trivial as exercising early or not will yield zero expected risk-neutral payoff. So exercising early isn't optimal in this case.

In all positions where the strike price is relative to $S,Su$, and $Sd$, we showed that exercising early is not optimal in all cases. Thus, early exercising an American-style call option (with no dividend yields) is never optimal.
