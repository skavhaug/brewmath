import numpy as np
import matplotlib.pylab as plt


def f(T):
    """
    :param T: Temperature in C
    :return: scale factor
    """
    return 2.39e11*np.exp(-9773/(273.15 + T))


def g(G):
    """
    :param G: Gravity of wort
    :return: bigness factor
    """
    return 1.65*1.25e-4**(G - 1)


def dUdt(t):
    return np.exp(-t/25)/(25*4.15)


def U(t):
    """
    :param t: boil time
    :return: boil time factor
    """
    return (1 - np.exp(-t/25))/4.15


def aaa(alpha_factor, hops, V):
    """
    :param alpha_factor: alpha acid content of hops (13% gives 0.13)
    :param hops: amount of hops in grams
    :param V: kettle volume in liters
    :return: mg/l of added alpha acid
    """
    return alpha_factor*hops*1000/V


def decimal_alpha_acid_util(G, t):
    """
    :param G: Gravity of wort
    :param t: Boil time
    :return: Utilization factor
    """
    return U(t)*g(G)


def tinseth(alpha_factor, hops, V, G, t):
    """
    :param alpha_factor: Alpha acid content of hops (13% gives 0.13)
    :param hops: Amount of hops in grams
    :param V: Kettle volume in liters
    :param G: Gravity of wort
    :param t: Boil time
    :return: IBUs
    """
    return decimal_alpha_acid_util(G, t)*aaa(alpha_factor, hops, V)


if __name__ == "__main__":
    alpha_factor = 0.13
    hops = 20
    V = 20
    G = 1.060
    t = np.linspace(0, 120, 121)
    plt.plot(t, tinseth(alpha_factor, hops, V, G, t))
    plt.show()
