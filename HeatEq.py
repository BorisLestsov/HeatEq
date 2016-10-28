import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors


class HeatEqSolver:
    def __init__(self, a=1.0, f=None,
                 x1=0.0, x2=1.0, t_fin=1.0,
                 phi=None, alpha=None, beta=None,
                 x_points=10, t_points=10,
                 method="expl", sigma = 0.75):

        assert x1 < x2
        assert t_fin > 0
        assert x_points > 0
        assert t_points > 0

        if method == "expl":
            self.solve = self._expl
        elif method == "impl":
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
        print self.tau, self.a ** 2, self.h ** 2, self.tau * (self.a ** 2) / (self.h ** 2)

    def _expl(self):
        self.sol = np.zeros(shape=(self.t_points, self.x_points))

        for i in xrange(self.x_points):
            self.sol[0, i] = self.phi(i*self.h)

        for t in xrange(1, self.t_points):
            self.sol[t, 0] = self.alpha(self.tau*t)
            self.sol[t, -1] = self.beta(self.tau*t)

        for t in xrange(1, self.t_points):
            for i in xrange(1, self.x_points-1):
                self.sol[t, i] = \
                    self.sol[t-1, i] + self.tau*(self.f(t-1, i*self.h) + (self.a**2)/(self.h**2) *
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
            self.sol[t, 0] = self.alpha(self.tau*t)
            self.sol[t, -1] = self.beta(self.tau*t)

        for t in xrange(1, self.t_points):
            self.m.fill(0)
            self.r.fill(0)
            self.m[0, 0] = 1
            self.r[0] = self.sol[t, 0]
            self.m[-1, -1] = 1
            self.r[-1] = self.sol[t, -1]
            self.r[0] = self.sol[t, 0]
            for i in xrange(1, self.x_points - 1):
                self.m[i, i - 1] = -self.utmp
                self.m[i, i] = 1 + 2*self.utmp
                self.m[i, i + 1] = -self.utmp
                self.r[i] = self.tau*self.f(t*self.tau, i*self.h) - self.dtmp*self.sol[t-1, i-1] + \
                            (1 + 2*self.dtmp)*self.sol[t-1, i] - self.dtmp*self.sol[t-1, i+1]

            #print self.m
            #print self.r
            #print '*'*60

            self.sol[t] = np.linalg.solve(self.m, self.r)



    def visualize(self, type="graph"):
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

        #cont = plt.pcolor((self.sol[t], self.sol[t]), cmap="rainbow")
        #cont = plt.plot(range(self.x_points), self.sol[t])

    def _empty(self, t = None, x = None):
        return 0


def main():
    a = 0.02
    x1 = 0.0
    x2 = 10.0
    x_points = 100
    t_fin = 8000.0
    t_points = 200

    def phi(x):
        return x/10 + np.math.exp(-(x-(x2-x1)/2)**2)

    def f(t, x):
        return 0.001*np.math.exp(-(x-(x2 - x1) / 5) ** 2) if t < t_fin/5 else 0

    def alpha(t):
        return 1

    def beta(t):
        return 0.8

    solver = HeatEqSolver(a=a, x1=x1, x2=x2, x_points=x_points,
                          t_fin=t_fin, t_points=t_points, phi=phi,
                          method="impl", f=f, alpha=alpha, beta=beta)

    solver.solve()
    solver.visualize(type="graph")


if __name__ == "__main__":
    main()