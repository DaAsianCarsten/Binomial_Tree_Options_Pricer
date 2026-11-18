import numpy as np 
import matplotlib.pyplot as plt 
from dataclasses import dataclass
from black_scholes_pricer import black_scholes_price

def call_payoff(K,S):
    return np.maximum(S-K,0.0)

def put_payoff(K,S):
    return np.maximum(K-S,0.0)

@dataclass
class binomial_option_pricer:
    # In years
    T: float
    # Number of Steps
    N: int
    S: float
    sigma: float
    r: float
    K: float
    call: bool=True
    european: bool=True
    q_d: float=0.0
    # Are we testing functions in the class?
    test: bool=False


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
        self.q = (np.exp((self.r - self.q_d) * self.dt) - self.d) / (self.u - self.d)

        # Check no-arbitrage rule
        if self.R <= self.d or self.R >= self.u:
            raise ValueError(f"No-Arbitrage Rule violated: {'R <= d' if self.R <= self.d else 'R >= u'}")

        # Initalized trees to be calculated
        self.option_price_tree = np.zeros((self.N+1, self.N+1))
        self.stock_price_tree = np.zeros((self.N+1, self.N+1))

        # Initialize our ending price
        self.V0 = None

        # For American options to show which nodes show be exercised early
        self.marked_nodes = np.zeros((self.N+1, self.N+1))

        # Generate the trees
        self.generate_trees()

        # Initialize and calculate Greeks
        self.delta = None
        self.gamma = None
        self.theta = None
        self.vega = None
        self.rho = None

        self.calculate_greeks()

        # Store parameters in dicts
        self.input_parameters = {"Time-span (yrs)": self.T, 
                                 "Time periods": self.N, 
                                 "Spot Price": self.S, 
                                 "Annualuzed Volatility": self.sigma, 
                                 "Continuous Risk-Free Rate": self.r,
                                 "Strike Price": self.K,
                                 "Option Type": ("European" if self.european else "American") + (" Call" if self.call else " Put"),
                                 "Dividend Yield": self.q_d}

        self.calculated_values = {"Time-step (yrs)": self.dt,
                                  "Risk-Free growth Factor": self.R,
                                  "Up-move Factor": self.u,
                                  "Down-move Factor": self.d,
                                  "Risk-Neutral Probability": self.q,
                                  "Fair Option Price": self.V0,
                                  "Delta (Δ)": self.delta,
                                  "Gamma (Γ)": self.gamma,
                                  "Theta (Θ)": self.theta,
                                  "Vega (v)": self.vega,
                                  "Rho (ρ)": self.rho}

    def verify_input_parameters(self):
        # All values should be finite
        if not all(map(np.isfinite,[self.T,self.N,self.S,self.sigma,self.r,self.K])): raise ValueError("All values must be finite")

        # T > 0
        if (self.T <= 0): raise ValueError("Time span T mist be positive")

        # N > 0 and integer
        if isinstance(self.N, (bool, np.bool_)) or not isinstance(self.N, (int, np.integer)): raise ValueError("# of time steps N must be a positive integer")

        # S > 0 
        if (self.S <= 0): raise ValueError("Spot price S must be positive")

        # K > 0
        if (self.K <= 0): raise ValueError("Strike price K must be positive")

        # sigma > 0
        if (self.sigma <= 0): raise ValueError("Annualized Volatility sigma must be positive")

        # Q_d >= 0
        if (self.q_d < 0): raise ValueError("Dividend Yield q_d must be nonnegative")


    def get_parameters(self,print_to_terminal=False):

        if print_to_terminal:
            print(f"""======================""")
            print(f"""Input Parameters: """)

            for k,v in self.input_parameters.items():
                print(f"{k}: {v}")

            print(f"""======================""")

            print(f"""Calculated Values: """)

            for k,v in self.calculated_values.items():
                print(f"{k}: {v}")

            print(f"""======================""")
        

        return self.input_parameters, self.calculated_values
        

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

    def calculate_greeks(self):
        # Delta
        V_u = self.option_price_tree[1,1]
        V_d = self.option_price_tree[1,0]
        S_u = self.stock_price_tree[1,1]
        S_d = self.stock_price_tree[1,0]
        self.delta = (V_u-V_d)/(S_u-S_d)

        # Gamma
        S_uu = self.stock_price_tree[2,2]
        S_ud = self.stock_price_tree[2,1]
        S_dd = self.stock_price_tree[2,0]
        V_uu = self.option_price_tree[2,2]
        V_ud = self.option_price_tree[2,1]
        V_dd = self.option_price_tree[2,0]

        S_u = self.stock_price_tree[1,1]
        S_d = self.stock_price_tree[1,0]

        Delta_u = (V_uu-V_ud)/(S_uu-S_ud)
        Delta_d = (V_ud-V_dd)/(S_ud-S_dd)
        self.gamma = (Delta_u-Delta_d)/(S_u-S_d)

        # Theta
        self.theta = (V_ud-self.V0)/(2*self.dt)

        # For the continuous calculations
        h = 1e-6
        try: assert(self.european or self.call)
        except:
            if not self.test:
                print("Unable to calculate Vega and Rho from continuous model due to American-Put option type.")
            return None
        
        # Vega (have to use closed form solution for this)
        V = lambda s: black_scholes_price(S=self.S,K=self.K,r=self.r,T=self.T,sigma=s,call=self.call,european=self.european)
        self.vega = (V(self.sigma+h)-V(self.sigma-h))/(2*h)

        # Rho (have to use closed form solution for this)
        V = lambda rr: black_scholes_price(S=self.S,K=self.K,r=rr,T=self.T,sigma=self.sigma,call=self.call,european=self.european)
        self.rho = (V(self.r+h)-V(self.r-h))/(2*h)




if __name__ == "__main__":
    S = 220
    K = 250
    T = 1
    N = 30
    sigma = 0.6
    r = 0.02
    call = True
    european = True
    q_d = 0.03

    A_with_dividend = binomial_option_pricer(T=T,N=N,S=S,sigma=sigma,r=r,call=call,european=european,q_d=q_d)

    A_with_dividend.get_parameters(print_to_terminal=True)