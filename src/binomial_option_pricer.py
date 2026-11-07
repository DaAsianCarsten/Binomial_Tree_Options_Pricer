import numpy as np 
import matplotlib.pyplot as plt 
from dataclasses import dataclass

def call_payoff(K,S):
    return np.maximum(S-K,0.0)

def put_payoff(K,S):
    return np.maximum(K-S,0.0)

@dataclass
class binomial_option_pricer:
    T: int=1
    # Number of Steps
    N: int=10
    S: float=100.0
    sigma: float=0.5
    r: float=0.02
    call: bool=True
    european: bool=True
    K: float=100


    def __post_init__(self):
        # Verify input parameters are safe to input
        self.verify_input_parameters()

        # Calculate the other parameters

        # Time step
        self.dt = self.T/self.N 

        # Risk-free growth factor
        self.R = np.exp(self.r*self.dt)

        # Up-move and down-move factor
        self.u = np.exp(self.sigma*np.sqrt(self.dt)) 
        self.d = 1/self.u

        # Risk-neutral probability
        self.q = (self.R - self.d) / (self.u - self.d)

        # Check no-arbitrage rule
        if self.R <= self.d or self.R >= self.u:
            raise ValueError(f"No-Arbitrage Rule violated: {'R <= d' if self.R <= self.d else 'R >= u'}")

        # Initalized trees to be calculated
        self.option_price_tree = np.zeros((self.N+1, self.N+1))
        self.stock_price_tree = np.zeros((self.N+1, self.N+1))

        # Initialize our ending price
        self.V0 = 0

        # For American options to show which nodes show be exercised early
        self.marked_nodes = np.zeros((self.N+1, self.N+1))

        # Generate the trees
        self.generate_trees()

        # Store parameters in dicts
        self.input_parameters = {"Time-span": self.T, 
                                 "Time periods": self.N, 
                                 "Spot Price": self.S, 
                                 "Annualuzed Volatility": self.sigma, 
                                 "Continuous Risk-Free Rate": self.r,
                                 "Strike Price": self.K,
                                 "Option Type": ("European" if self.european else "American")+ (" Call" if self.call else " Put")}

        self.calculated_parameters = {"Time-step": self.dt,
                                      "Risk-Free growth Factor": self.R,
                                      "Up-move Factor": self.u,
                                      "Down-move Factor": self.d,
                                      "Risk-Neutral Probability": self.q,
                                      "Fair Option Price": self.V0}

    def verify_input_parameters(self):
        # All values should be finite
        if not all(map(np.isfinite,[self.T,self.N,self.S,self.sigma,self.r,self.K])): raise ValueError("All values must be finite")

        # T > 0
        if (self.T <= 0): raise ValueError("Time span T mist be positive")

        # N > 0 and integer
        if (self.N <= 0 or self.N != int(self.N)): raise ValueError("# of time steps N must be a positive integer")

        # S > 0 
        if (self.S <= 0): raise ValueError("Spot price S must be positive")

        # K > 0
        if (self.K <= 0): raise ValueError("Strike price K must be positive")

        # sigma > 0
        if (self.sigma <= 0): raise ValueError("Annualized Volatility sigma must be positive")


    def get_parameters(self):
        print(
        f"""
        ===============================
        Input Parameters:
        Time-span (years) T = {self.T}
        Time periods N = {self.N}
        Spot Price S = {self.S}
        Annualized Volatility sigma = {self.sigma} 
        Continuous Risk-Free Gain r = {self.r}
        Strike Price K = {self.K}
        Option Type: {("European" if self.european else "American")+ (" Call" if self.call else " Put")}
        ===============================
        Calculated Parameters:
        Time-step dt = {self.dt} 
        Continuous Risk-Free Growth Factor R = {self.R} 
        Up-Move Gain u = {self.u} 
        Down-Move Loss d = {self.d}
        Risk-Neutral Probability q = {self.q}
        Fair Option Price V0 = {round(self.V0,2)}
        ===============================
        """)

        return self.input_parameters, self.calculated_parameters
        

    def generate_trees(self):
        g = call_payoff if self.call else put_payoff
        V = g(self.K, np.array([self.S * self.u**j * self.d**(self.N-j) for j in range(self.N+1)]))
        self.option_price_tree[self.N] = V
        self.stock_price_tree[self.N] = np.array([self.S * self.u**j * self.d**(self.N-j) for j in range(self.N+1)])

        for i in range(self.N-1,-1,-1):
            V[:i+1] = 1/self.R * (self.q * V[1:i+2] + (1-self.q) * V[:i+1])
            if not self.european: 
                b_array = (V[:i+1] < g(self.K, np.array([self.S * self.u**j * self.d**(i-j) for j in range(i+1)])))
                V[:i+1] = np.maximum(V[:i+1], g(self.K, np.array([self.S * self.u**j * self.d**(i-j) for j in range(i+1)])))

                self.marked_nodes[i,np.where(b_array)[0]] = 1
            V[i+1:] = 0

            self.option_price_tree[i] = V
            self.stock_price_tree[i,:i+1] = np.array([self.S * self.u**j * self.d**(i-j) for j in range(i+1)])

        self.V0 = self.option_price_tree[0][0]



if __name__ == "__main__":
    pass