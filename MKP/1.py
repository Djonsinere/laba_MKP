import math
import numpy as np
import matplotlib.pyplot as plt


e = 0.0549   # Эксцентриситет
T = 27.321661*24*3600 # Период обращения в секундах (один оборот)
a = 384748   # Большая полуось 
t = np.linspace(0, T, 1000) # Временные интервалы для одного оборота (0 до T)
M = t / T * 2* np.pi

#Метод итераций (метод последовательных приближений)
def equation_iterations(M, e, tol=1e-6):
    E = M
    E_new = e * np.sin(E) + M
    while E_new - E > tol:
        E = E_new
        E_new = e * np.sin(E) + M
    return E_new

#Метод золотого сечения
def golden_section(M, e, tol=1e-6):
    def f(E):
        return E - e * np.sin(E) - M
    a, b = 0, 2 * np.pi
    phi = (1 + np.sqrt(5)) / 2
    while (b - a) > tol:
        c = a + (b - a) / phi
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

#Функция для вычисления эксцентрической аномалии с использованием метода Ньютона
def equation_newton(M, e):
    E = M
    for _ in range(100):
        E_new = E + (M - (E - e * np.sin(E))) / (1 - e * np.cos(E))
        if abs(E_new - E) < 1e-6:
            break
        E = E_new
    return E



#Метод половинного деления
def equation_bisection(M, e, tol=1e-6):
    def f(E):
        return E - e * np.sin(E) - M
    a, b = 0, 2 * np.pi
    while (b - a) > tol:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

methods = {
    "Ньютона": equation_newton,
    "Золотого сечения": golden_section,
    "Половинного деления": equation_bisection,
    "Итераций": equation_iterations,
}

results = {}
timing_results = {}

for name, method in methods.items():
    E = np.array([method(m, e) for m in M])
    results[name] = E
    
# Вывод результатов
for name, E in results.items():
    print(f"Метод {name}: E = {E[-1]:.8f}")

#Функция для вычисления истинной аномалии
def true_anomaly(E, e):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

# Вычисление истинной аномалии для метода Ньютона
nu = np.array([true_anomaly(equation_newton(m, e), e) for m in M])  
E_newton = results["Ньютона"]
# Коррекция аномалий, чтобы избежать падения в ноль в конце
if M[-1] < 1e-6:
    M[-1] = M[-2]
if E_newton[-1] < 1e-6:
    E_newton[-1] = E_newton[-2]
if nu[-1] < 1e-6:
    nu[-1] = nu[-2]

# Построение графиков
plt.figure(figsize=(12, 8))
plt.plot(t, M, label='Средняя аномалия M(t)', color='blue')
plt.plot(t, E_newton, label='Эксцентрическая аномалия E(t)', color='red')
plt.plot(t, nu, label='Истинная аномалия ν(t)', color='green')
plt.title('Зависимости аномалий от времени для 1 оборота Lunar Reconnaissance Orbiter (метод Ньютона)')
plt.xlabel('Время')
plt.ylabel('Аномалия')
plt.legend()
plt.grid()
plt.tight_layout()
plt.yticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi], ['0', 'π/2', 'π', '3π/2', '2π'])
plt.show()

#Кусочки кода для построения графиков другими методами
def extra_code(e, M):
    # Вычисление истинной аномалии для метода золотого сечения
    gs = np.array([true_anomaly(golden_section(m, e), e) for m in M])
    E_newton = results["Золотого сечения"]
    # Коррекция аномалий, чтобы избежать падения в ноль в начале
    if M[0] < 1e-6 or E_newton[0] < 1e-6 or gs[0] < 1e-6:
        M[0] = M[1]
        E_newton[0] = E_newton[1]
        gs[0] = gs[1]
    # Вычисление истинной аномалии для метода итераций
    it = np.array([true_anomaly(equation_iterations(m, e), e) for m in M])
    E_newton = results["Итераций"]
    # Коррекция аномалий, чтобы избежать падения в ноль в конце
    if M[-1] < 1e-6:
        M[-1] = M[-2]
    if E_newton[-1] < 1e-6:
        E_newton[-1] = E_newton[-2]
    if it[-1] < 1e-6:
        it[-1] = it[-2]
    # Вычисление истинной аномалии для метода половинного деления
    bs = np.array([true_anomaly(equation_bisection(m, e), e) for m in M])
    E_newton = results["Половинного деления"]
    # Коррекция аномалий, чтобы избежать падения в ноль в начале
    if M[0] < 1e-6 or E_newton[0] < 1e-6 or bs[0] < 1e-6:
        M[0] = M[1]
        E_newton[0] = E_newton[1]
        bs[0] = bs[1]