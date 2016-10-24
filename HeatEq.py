import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class HeatEqSolver:
    def __init__(self, a=1, f=None,
                 x1=0, x2=1, t_fin=1,
                 alpha1=None, alpha2=None,
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
        if alpha1 is not None:
            self.alpha1 = alpha1
        else:
            self.alpha1 = self._empty
        if alpha1 is not None:
            self.alpha2 = alpha2
        else:
            self.alpha2 = self._empty
        self.x_points = x_points
        self.t_points = t_points
        self.h = float(x2 - x1)/x_points
        self.tau = float(t_fin) / t_points
        self.arr = np.zeros(shape=(t_points, x_points))
        self.arr = np.random.randn(self.x_points, self.t_points)
    def solve(self):
        pass

    def visualize(self):
        def init():
            x = range(self.x_points)
            plt.xlabel(x)
            plt.ylabel(None)
            cont = plt.contourf((self.arr[0], self.arr[0]))
            return cont,

        def animate(t):
            cont = plt.contourf((self.arr[t], self.arr[t]))
            return cont,

        fig = plt.figure(figsize=(80, 5), dpi=10)
        ani = animation.FuncAnimation(fig, animate, frames=self.t_points, interval=1,
                                      repeat=False, init_func=init, )
        plt.show()


    def _empty(self, t = None, x = None):
        return 0

def main():
    solver = HeatEqSolver(x_points=50, t_points=50)

    solver.solve()
    solver.visualize()





if __name__ == "__main__":
    main()