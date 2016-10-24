import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class HeatEqSolver:
    def __init__(self, a=1, f=None,
                 x1=0, x2=1, t_fin=1,
                 phi=None, alpha=None, beta=None,
                 x_points=10, t_points=10):

        assert x1 < x2
        assert t_fin > 0
        assert x_points > 0
        assert t_points > 0

        self.a = float(a)
        if f is not None:
            self.f = f
        else:
            self.f = self._empty
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.t_fin = float(t_fin)
        if alpha is not None:
            self.alpha = alpha
        else:
            self.alpha = self._noalpha
        if beta is not None:
            self.beta = beta
        else:
            self.beta = self._nobeta
        if phi is not None:
            self.phi = phi
        else:
            self.phi = self._empty
        self.x_points = x_points
        self.t_points = t_points
        self.h = float(x2 - x1)/x_points
        self.tau = float(t_fin) / t_points

    def solve(self):
        self.sol = np.zeros(shape=(self.t_points, self.x_points))

        for i in xrange(1, self.x_points):
            self.sol[0, i] = self.phi(i)

        for t in xrange(1, self.t_points):
            self.sol[t, 0] = self.alpha(t)
            self.sol[t, self.x_points-1] = self.beta(t)
            for i in xrange(1, self.x_points-1):
                self.sol[t, i] = \
                    self.sol[t-1, i] + self.tau * (self.f(t-1, i) + (self.a**2)/(self.h**2) *
                    (self.sol[t-1, i-1] - 2*self.sol[t-1, i] + self.sol[t-1, i+1] + self.f(t-1, i)))

    def visualize(self):
        def init():
            x = range(self.x_points)
            plt.xlabel(x)
            plt.ylabel(None)
            return plt

        def animate(t):
            cont = plt.pcolor((self.sol[t], self.sol[t]))
            return cont

        fig = plt.figure(figsize=(80, 5), dpi=10)
        ani = animation.FuncAnimation(fig, animate, frames=self.t_points, interval=200,
                                      repeat=False, init_func=init)
        plt.show()

    def _empty(self, t = None, x = None):
        return 0

    def _noalpha(self, t):
        return self.sol[t, 1]

    def _nobeta(self, t):
        return self.sol[t, -2]


def main():
    a = 0.015
    x1 = 0
    x2 = 1
    x_points = 500
    t_fin = 1
    t_points = 100

    def phi(i):
        return i
    #def f(t, i):
    #    return x_points/2-abs(x_points/2-i)
    def ff(x):
        return 0

    solver = HeatEqSolver(a=a, x1=x1, x2 = x2, x_points=x_points,
                          t_fin=t_fin, t_points=t_points, phi=phi, alpha=ff, beta=ff)

    solver.solve()
    solver.visualize()



if __name__ == "__main__":
    main()