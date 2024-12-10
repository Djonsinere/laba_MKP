import numpy as np
import matplotlib.pyplot as plt

T = 27.321661*24*3600  # Период обращения в секундах (один оборот)
a = 384748   # Большая полуось в а.е.
mu = 4902.800066  # Гравитационный параметр Луны в км^3/с^2
e = 0.0549  # Эксцентриситет

#Метод Ньютона для решения уравнения Кеплера
def Newton(M, e):
    E = M
    for _ in range(100):
        E_new = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        if abs(E_new - E) < 1e-6:
            break
        E = E_new
    return E

#Расчет средней аномалии M(t)
def calculate_mean_anomaly(t, T):
    n = 2 * np.pi / T
    M = n * t
    return M

#Расчет истинной аномалии ν(t)
def true_anomaly(E, e):
    beta = np.sqrt((1 + e) / (1 - e))
    nu = 2 * np.arctan(beta * np.tan(E / 2))
    return nu

#Расчет эксцентрической аномалии E(t)
def eccentric_anomaly(M, e):
    E = np.array([Newton(m, e) for m in M])
    return E

#Расчет радиальной скорости Vr(t)
def radial_velocity(nu, a, e, mu):
    p = a * (1 - e**2)
    Vr = np.sqrt(mu / p) * e * np.sin(nu)
    return Vr

#Расчет радиус-вектора r(t)
def radius_vector(nu, a, e):
    p = a * (1 - e**2)
    r = p / (1 + e * np.cos(nu))
    return r

#модуля скорости V(t)
def calculate_speed(Vn, Vr):
    V = np.sqrt(Vn**2 + Vr**2)
    return V
    
#Расчет трансверсальной скорости Vn(t)
def transversal_velocity(nu, a, e, mu):
    p = a * (1 - e**2)
    Vn = np.sqrt(mu / p) * (1 + e * np.cos(nu))
    return Vn

def plot_graphs(t, r, Vr, Vn, V):
    plt.figure(figsize=(10, 10))
    # радиус-вектор
    plt.subplot(2, 2, 1)
    plt.plot(t, r, label='r (радиус-вектор)', color='blue')
    plt.title('Радиус-вектор')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('r (м)')
    plt.grid()
    plt.legend()
    # радиальной скорость
    plt.subplot(2, 2, 2)
    plt.plot(t, Vr, label='Vr(t) (радиальная скорость)', color='red')
    plt.title('Радиальная скорость')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('Vr (м/с)')
    plt.grid()
    plt.legend()
    # трансверсальной скорость
    plt.subplot(2, 2, 3)
    plt.plot(t, Vn, label='Vn(t) (трансверсальная скорость)', color='orange')
    plt.title('Трансверсальная скорость')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('Vn (м/с)')
    plt.grid()
    plt.legend()
    # модуля скорость
    plt.subplot(2, 2, 4)
    plt.plot(t, V, label='V (модуль скорости)', color='pink')
    plt.title('Модуль скорости')
    plt.xlabel('Время t (в секундах)')
    plt.ylabel('V (м/с)')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()

def main():
    t = np.linspace(0, T, 1000)
    M = calculate_mean_anomaly(t, T)
    E = eccentric_anomaly(M, e)
    nu = true_anomaly(E, e)
    r = radius_vector(nu, a, e)
    Vr = radial_velocity(nu, a, e, mu)
    Vn = transversal_velocity(nu, a, e, mu)
    V = calculate_speed(Vn, Vr)
    plot_graphs(t / (24 * 3600), r, Vr, Vn, V)  #в сутки

if __name__ == "__main__":
    main()

