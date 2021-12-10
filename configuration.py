""" configuration object """


class OptionConfigurationBuilder(object):
    """ configuration object """
    def __init__(self, kind: str = None, spot: float = None, strike: float = None, sigma: float = None,
                 maturity: int = None, risk_free_rate: float = None, dividend_yield: float = None,
                 simulation: int = None, steps: int = None, barrier: float = None):
        """
        :param kind: should be 'call' or 'put'
        :param spot: current price
        :param strike: exercise price
        :param sigma: current volatility
        :param maturity: number of day until expiration
        :param risk_free_rate: theoretical rate of return of an investment with zero risk
        :param dividend_yield: financial ratio (dividend/price)
        :param simulation: number of montecarlo simulation
        :param steps: number of steps in each simulation
        :param barrier: barrier for exotic options
        """
        self._kind = kind
        self._spot = spot
        self._strike = strike
        self._sigma = sigma
        self._maturity = maturity
        self._risk_free_rate = risk_free_rate
        self._dividend_yield = dividend_yield
        self._simulation = simulation
        self._steps = steps
        self._barrier = barrier

    @property
    def kind(self):
        """ getter """
        if self._kind not in ['call', 'put']:
            raise ValueError('kind should be call or put')
        return self._kind

    @kind.setter
    def kind(self, _kind):
        """ setter """
        if _kind not in ['call', 'put']:
            raise ValueError('kind should be call or put')
        self._kind = _kind

    @property
    def spot(self):
        """ getter """
        if not isinstance(self._spot, (float, int)):
            raise ValueError('spot should be a float or an integer')
        return self._spot

    @spot.setter
    def spot(self, _spot):
        """ setter """
        if not isinstance(_spot, (float, int)):
            raise ValueError('spot should be a float or an integer')
        self._spot = _spot

    @property
    def strike(self):
        """ getter """
        if not isinstance(self._strike, (float, int)):
            raise ValueError('strike should be a float or an integer')
        return self._strike

    @strike.setter
    def strike(self, _strike):
        """ setter """
        if not isinstance(_strike, (float, int)):
            raise ValueError('strike should be a float or an integer')
        self._strike = _strike

    @property
    def sigma(self):
        """ getter """
        if not isinstance(self._sigma, (float, int)):
            raise ValueError('sigma should be a float or an integer')
        return self._sigma

    @sigma.setter
    def sigma(self, _sigma):
        """ setter """
        if not isinstance(_sigma, (float, int)):
            raise ValueError('sigma should be a float or an integer')
        self._sigma = _sigma

    @property
    def maturity(self):
        """ getter """
        if not isinstance(self._maturity, int):
            raise ValueError('maturity should be an integer')
        return self._maturity

    @maturity.setter
    def maturity(self, _maturity):
        """ setter """
        if not isinstance(_maturity, int):
            raise ValueError('maturity should be an integer')
        self._maturity = _maturity

    @property
    def risk_free_rate(self):
        """ getter """
        if not isinstance(self._risk_free_rate, (int, float)):
            raise ValueError('risk_free_rate should be a float or an integer')
        return self._risk_free_rate

    @risk_free_rate.setter
    def risk_free_rate(self, _risk_free_rate):
        """ setter """
        if not isinstance(_risk_free_rate, (int, float)):
            raise ValueError('risk_free_rate should be a float or an integer')
        self._risk_free_rate = _risk_free_rate

    @property
    def dividend_yield(self):
        """ getter """
        if not isinstance(self._dividend_yield, (int, float)):
            raise ValueError('dividend_yield should be a float or an integer')
        return self._dividend_yield

    @dividend_yield.setter
    def dividend_yield(self, _dividend_yield):
        """ setter """
        if not isinstance(_dividend_yield, (int, float)):
            raise ValueError('dividend_yield should be a float or an integer')
        self._dividend_yield = _dividend_yield

    @property
    def simulation(self):
        """ getter """
        if not isinstance(self._simulation, int):
            raise ValueError('simulation should be an integer')
        return self._simulation

    @simulation.setter
    def simulation(self, _simulation):
        """ setter """
        if not isinstance(_simulation, int):
            raise ValueError('simulation should be an integer')
        self._simulation = _simulation

    @property
    def steps(self):
        """ getter """
        if not isinstance(self._steps, int):
            raise ValueError('steps should be an integer')
        return self._steps

    @steps.setter
    def steps(self, _steps):
        """ setter """
        if not isinstance(_steps, int):
            raise ValueError('steps should be an integer')
        self._steps = _steps

    @property
    def barrier(self):
        """ getter """
        if self._barrier is not None and self._barrier is not isinstance(self._barrier, (int, float)):
            raise ValueError('barrier should be a float or an integer')
        return self._barrier

    @barrier.setter
    def barrier(self, _barrier):
        """ setter """
        if not isinstance(_barrier, int):
            raise ValueError('barrier should be a float or an integer')
        self._barrier = _barrier
