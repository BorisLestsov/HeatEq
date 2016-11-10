import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors


class HeatEqSolver:
    def __init__(self, a=1.0, f=None,
                 x1=0.0, x2=1.0, t_fin=1.0,
                 phi=None, alpha=None, beta=None,
                 x_points=10, t_points=10,
                 ab_cond="first",
                 bb_cond="first",
                 method="explicit", sigma = 0.75):
        """
        This method solves one dimensional Heat Equation
        du/dt = a^2*d2u/d2t + f
        with boundary conditions:
        u(0, x) = phi(x)
        u(t, 0) = alpha(t)
        u(t, M) = beta(t)

        :param a:           Heat conductivity coefficient
        :param f:           f fucntion in the equation
        :param x1:          left bound
        :param x2:          right bound
        :param t_fin:       time bound
        :param phi:         initial conditions
        :param alpha:       left boundary condition
        :param beta:        right boundary condition
        :param x_points:    number of points in grid for x
        :param t_points:    number of points in grid for t
        :param b_cond:      type of boundary conditions ("first" or "second")
        :param method:      explicit of implicit solver type ("explicit" of "implicit")
        :param sigma:       sigma parameter in implicit solver
        """
        assert x1 < x2
        assert t_fin > 0
        assert x_points > 0
        assert t_points > 0

        if method == "explicit":
            self.solve = self._expl
        elif method == "implicit":
            self.solve = self._impl
            assert sigma <= 1 and sigma >= 0
            self.sigma = sigma
        else:
            raise Exception("Wrong method")

        self.a = a
        if f is not None:
            self.f = f
        else:
            self.f = self._empty
        self.x1 = x1
        self.x2 = x2
        self.t_fin = t_fin
        if alpha is not None:
            self.alpha = alpha
        else:
            self.alpha = self._empty
        if beta is not None:
            self.beta = beta
        else:
            self.beta = self._empty
        if phi is not None:
            self.phi = phi
        else:
            self.phi = self._empty
        self.x_points = x_points
        self.t_points = t_points
        self.h = (x2 - x1)/x_points
        self.tau = t_fin / t_points
        if ab_cond == "first":
            self.ab_cond = "first"
        elif ab_cond == "second":
            self.ab_cond = "second"
        else:
            raise Exception("Wrong boundary condition")
        if bb_cond == "first":
            self.bb_cond = "first"
        elif bb_cond == "second":
            self.bb_cond = "second"
        else:
            raise Exception("Wrong boundary condition")
        print self.tau, self.a ** 2, self.h ** 2, self.tau * (self.a ** 2) / (self.h ** 2)

    def _expl(self):
        self.sol = np.zeros(shape=(self.t_points, self.x_points))

        for i in xrange(self.x_points):
            self.sol[0, i] = self.phi(i*self.h)

        #for t in xrange(1, self.t_points):
        #    self.sol[t, 0] = self.alpha(self.tau*t)
        #   self.sol[t, -1] = self.beta(self.tau*t)

        for t in xrange(1, self.t_points):
            if self.ab_cond == "first":
                self.sol[t, 0] = self.alpha(self.tau * t)
            else:
                self.sol[t, 0] = -self.alpha(self.tau * t) * self.h + self.sol[t, 1]
            if self.bb_cond == "first":
                self.sol[t, -1] = self.beta(self.tau * t)
            else:
                self.sol[t, -1] = self.alpha(self.tau * t) * self.h + self.sol[t, -2]
            for i in xrange(1, self.x_points-1):
                self.sol[t, i] = \
                    self.sol[t-1, i] + self.tau*(self.f(t*self.tau, i*self.h) + (self.a**2)/(self.h**2) *
                    (self.sol[t-1, i-1] - 2*self.sol[t-1, i] + self.sol[t-1, i+1]))


    def _impl(self):
        self.m = np.zeros(shape=(self.x_points, self.x_points))
        self.r = np.zeros(shape=(self.x_points,))
        self.sol = np.zeros(shape=(self.t_points, self.x_points))

        self.utmp = (self.sigma * self.tau * self.a ** 2) / (self.h ** 2)
        self.dtmp = ((1 - self.sigma) * self.tau * self.a ** 2) / (self.h ** 2)

        for i in xrange(self.x_points):
            self.sol[0, i] = self.phi(i*self.h)

        for t in xrange(1, self.t_points):
            self.m.fill(0)
            self.r.fill(0)
            self.m[0, 0] = 1
            self.r[0] = self.sol[t, 0]
            self.m[-1, -1] = 1
            self.r[-1] = self.sol[t, -1]
            self.r[0] = self.sol[t, 0]
            if self.b_cond == "first":
                self.sol[t, 0] = self.alpha(self.tau * t)
                self.sol[t, -1] = self.beta(self.tau * t)
            for i in xrange(1, self.x_points - 1):
                self.m[i, i - 1] = -self.utmp + self.alpha(self.tau * t) * self.h + self.sol[t, 1]
                self.m[i, i] = 1 + 2*self.utmp + self.alpha(self.tau * t) * self.h + self.sol[t, -2]
                self.m[i, i + 1] = -self.utmp
                self.r[i] = self.tau*self.f(t*self.tau, i*self.h) - self.dtmp*self.sol[t-1, i-1] + \
                            (1 + 2*self.dtmp)*self.sol[t-1, i] - self.dtmp*self.sol[t-1, i+1]


            self.sol[t] = np.linalg.solve(self.m, self.r)



    def visualize(self, type="graph"):
        """
        Plot graph
        :param type: "graph" or "pcolor" - type of visualization
        :return: None
        """

        def _graph_animate(t):
            line.set_ydata(self.sol[t])
            return line,

        def _pcolor_animate(t):
            cont = plt.pcolor((self.sol[t], self.sol[t]),
                              norm=colors.Normalize(vmin=np.min(self.sol), vmax=np.max(self.sol)),
                              cmap="plasma")
            return cont

        x = [self.x1 + dx*self.h for dx in range(self.x_points)]
        if type == "graph":
            fig, ax = plt.subplots()
            line, = ax.plot(x, self.sol[0])
            ani = animation.FuncAnimation(fig, _graph_animate, frames=self.t_points, interval=5000.0/self.t_points, repeat=False)
        elif type == "pcolor":
            fig = plt.figure(figsize=(80, 5), dpi=10)
            ani = animation.FuncAnimation(fig, _pcolor_animate, frames=self.t_points, interval=10.0/self.t_points*100, repeat=False)
        else:
            raise Exception("Unknown plot type")

        plt.show()

    def _empty(self, t = None, x = None):
        return 0


def main():
    a = 0.01
    x1 = 0.0
    x2 = 1.0
    x_points = 100
    t_fin = 80.0
    t_points = 1000

    def phi(x):
        return np.math.exp(-(x-(x2-x1)/2)**2)


    def f(t, x):
        return 0.001*np.math.exp(-(x-(x2 - x1) / 5) ** 2) if t < t_fin/10 else -0.00002

    def alpha(t):
        return 1

    def beta(t):
        return 0.06*np.math.sin(t)

    solver = HeatEqSolver(a=a, x1=x1, x2=x2, x_points=x_points,
                          t_fin=t_fin, t_points=t_points, phi=None,
                          method="explicit", ab_cond="second", f=None,
                          alpha=alpha, beta=beta, bb_cond="first")

    solver.solve()
    solver.visualize(type="graph")


if __name__ == "__main__":
    main()