""" pricing module """
import numpy as np
import scipy.stats as stats
from math import log, sqrt, exp
from configuration import OptionConfigurationBuilder


class Option(object):
    """ Option object which will be use for addition, subtraction and multiplication """
    def __init__(self, price: float = None, delta: float = None, gamma: float = None, vega: float = None,
                 theta: float = None, rho: float = None):
        self._price = price
        self._delta = delta
        self._gamma = gamma
        self._vega = vega
        self._theta = theta
        self._rho = rho

    def price(self):
        """ Premium, composed of the sum of the intrinsic and time value """
        return self._price

    def delta(self):
        """ first derivative of the value of the option with respect to the underlying security's price """
        return self._delta

    def gamma(self):
        """ second derivative of the value of the option with respect to the underlying security's price """
        return self._gamma

    def vega(self):
        """ first derivative of the value of the option with respect to the underlying security's volatility """
        return self._vega

    def theta(self):
        """ first derivative of the value of the option with respect to the time """
        return self._theta

    def rho(self):
        """ first derivative of the value of the option with respect to the interest rate """
        return self._rho

    def __add__(self, other):

        return Option(
            self.price() + other.price(),
            self.delta() + other.delta(),
            self.gamma() + other.gamma(),
            self.vega() + other.vega(),
            self.theta() + other.theta(),
            self.rho() + other.rho()
        )

    def __sub__(self, other):

        return Option(
            self.price() - other.price(),
            self.delta() - other.delta(),
            self.gamma() - other.gamma(),
            self.vega() - other.vega(),
            self.theta() - other.theta(),
            self.rho() - other.rho()
        )

    def __mul__(self, other):

        return Option(
            self.price() * other,
            self.delta() * other,
            self.gamma() * other,
            self.vega() * other,
            self.theta() * other,
            self.rho() * other
        )

    def __hash__(self):
        return hash((self.price(), self.gamma(), self.vega(), self.theta(), self.rho()))


class BlackScholesMerton(Option):
    """ Black Scholes Merton Pricing """
    PERIODS_PER_YEAR = 252

    def __init__(self, configuration: OptionConfigurationBuilder):
        super(BlackScholesMerton, self).__init__()
        self.kind = configuration.kind
        self.s = configuration.spot
        self.k = configuration.strike
        self.v = configuration.sigma
        self.t = configuration.maturity / self.PERIODS_PER_YEAR
        self.r = configuration.risk_free_rate
        self.q = configuration.dividend_yield
        self._d1 = (log(self.s / self.k) + (self.r - self.q + self.v ** 2 * 0.5) * self.t) / (self.v * sqrt(self.t))
        self._d2 = (log(self.s / self.k) + (self.r - self.q - self.v ** 2 * 0.5) * self.t) / (self.v * sqrt(self.t))

    def price(self):
        """ Premium, composed of the sum of the intrinsic and time value """
        if self.kind == 'call':
            price = exp(-self.r * self.t) * (self.s * exp((self.r - self.q) * self.t) * stats.norm.cdf(
                self._d1) - self.k * stats.norm.cdf(self._d2))
        else:
            price = exp(-self.r * self.t) * (self.k * stats.norm.cdf(-self._d2) - (
                    self.s * exp((self.r - self.q) * self.t) * stats.norm.cdf(-self._d1)))
        return price

    def delta(self):
        """ first derivative of the value of the option with respect to the underlying security's price """
        if self.kind == 'call':
            delta = exp(-self.q * self.t) * stats.norm.cdf(self._d1)
        else:
            delta = exp(-self.q * self.t) * stats.norm.cdf(self._d1) - 1
        return delta

    def gamma(self):
        """ second derivative of the value of the option with respect to the underlying security's price """
        return stats.norm.pdf(self._d1) * exp(-self.q * self.t) / (self.s * self.v * sqrt(self.t))

    def vega(self):
        """ first derivative of the value of the option with respect to the underlying security's volatility """
        return 0.01 * (self.s * sqrt(self.t) * stats.norm.pdf(self._d1) * exp(-self.q * self.t))

    def theta(self):
        """ first derivative of the value of the option with respect to the time """
        if self.kind == 'call':
            theta = -self.s * stats.norm.pdf(self._d1) * self.v * exp(-self.q * self.t) / (2 * sqrt(self.t)) \
                    + self.q * self.s * stats.norm.cdf(self._d1) * exp(-self.q * self.t) \
                    - self.r * self.k * exp(-self.r * self.t) * stats.norm.cdf(self._d2)
        else:
            theta = -self.s * stats.norm.pdf(self._d1) * self.v * exp(-self.q * self.t) / (2 * sqrt(self.t)) \
                    + self.q * self.s * stats.norm.cdf(-self._d1) * exp(-self.q * self.t) \
                    - self.r * self.k * exp(-self.r * self.t) * stats.norm.cdf(-self._d2)
        return 1/self.PERIODS_PER_YEAR * theta

    def rho(self):
        """ first derivative of the value of the option with respect to the interest rate """
        if self.kind == 'call':
            rho = 0.01 * (self.k * self.t * (exp(-self.r * self.t)) * stats.norm.cdf(self._d2))
        else:
            rho = 0.01 * (-self.k * self.t * (exp(-self.r * self.t)) * stats.norm.cdf(-self._d2))
        return rho


class GeometricBrownianMotion(BlackScholesMerton):
    """ Continuous-time stochastic process """
    def __init__(self, configuration: OptionConfigurationBuilder):
        super(GeometricBrownianMotion, self).__init__(configuration)
        self.simulation = configuration.simulation
        self.steps = configuration.steps
        self.st_paths = np.zeros((self.simulation, self.steps))
        self.prices_at_maturity = None
        self.barrier = configuration.barrier

        GeometricBrownianMotion.run_simulation(self)

    def run_simulation(self):
        """ time series computation """
        for i in range(self.simulation):
            self.st_paths[i][0] = self.s
            for j in range(1, self.steps):
                self.st_paths[i][j] = self.st_paths[i][j-1] * exp(
                    (self.r - 0.5 * self.v**2) * 1 / self.PERIODS_PER_YEAR + self.v * sqrt(
                        1 / self.PERIODS_PER_YEAR) * np.random.normal(0, 1))
        self.prices_at_maturity = [self.st_paths[i][-1] for i in range(self.simulation)]

    def price(self):
        """ average of discounted payoffs """
        if self.kind == 'call':
            payoffs = [max(S - self.k, 0) for S in self.prices_at_maturity]
        else:
            payoffs = [max(self.k - S, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * exp(-self.r * self.t)

    def digital(self):
        """ average of discounted payoffs """
        payoffs = [S > self.k for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * exp(-self.r * self.t)

    def up_and_out(self):
        """ average of discounted payoffs """
        payoffs = [(S < self.barrier) * max(S - self.k, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)

    def up_and_in(self, barrier):
        """ average of discounted payoffs """
        payoffs = [(S > barrier) * max(S - self.k, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * exp(-self.r * self.t)

    def down_and_out(self, barrier):
        """ average of discounted payoffs """
        payoffs = [(S > barrier) * max(self.k - S, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * exp(-self.r * self.t)

    def down_and_in(self, barrier):
        """ average of discounted payoffs """
        payoffs = [(S < barrier) * max(self.k - S, 0) for S in self.prices_at_maturity]
        payoff = np.mean(payoffs)
        return payoff * exp(-self.r * self.t)