"""
Лабораторная работа № 1. Вариант 19
Итерационные методы вычисления sqrt(a), cbrt(b) и ln(c).

Выполнила: Ситникова Ольга
Группа: БИ2-1
Параметры варианта: a=53, b=300, c=6, eps=1e-4.
Для квадратного и кубического корней используется критерий по невязке.
Для натурального логарифма остановка выполняется по модулю очередного члена ряда.
Дополнительное требование Д1: вывод таблиц итераций.
"""

import math

A = 53.0
B = 300.0
C = 6.0
EPS = 1e-4
N_MAX = 10_000


def sqrt_newton(a: float, eps: float = EPS, n_max: int = N_MAX):
    """Вычисление sqrt(a) методом Ньютона; остановка по |x_new**2-a| < eps."""
    if a < 0:
        raise ValueError("Для квадратного корня аргумент a должен быть неотрицательным.")
    if eps <= 0:
        raise ValueError("Точность eps должна быть положительной.")
    if a == 0:
        return 0.0, 0, []
    x = float(a)  # По варианту 19: x0 = a
    trace = []
    for n in range(1, n_max + 1):
        x_new = 0.5 * (x + a / x)
        delta = abs(x_new - x)
        residual = abs(x_new * x_new - a)
        trace.append({
            "n": n, "x_n": x, "x_next": x_new,
            "delta": delta, "residual": residual
        })
        if residual < eps:
            return x_new, n, trace
        x = x_new
    raise RuntimeError(f"Точность не достигнута за {n_max} итераций.")


def cbrt_newton(b: float, eps: float = EPS, n_max: int = N_MAX):
    """Вычисление действительного кубического корня методом Ньютона."""
    if eps <= 0:
        raise ValueError("Точность eps должна быть положительной.")
    if b == 0:
        return 0.0, 0, []
    sign = -1.0 if b < 0 else 1.0
    value = abs(float(b))
    x = value  # Для варианта 19 b положительно; x0 = |b|
    trace = []
    for n in range(1, n_max + 1):
        x_new = (2.0 * x + value / (x * x)) / 3.0
        delta = abs(x_new - x)
        residual = abs(x_new ** 3 - value)
        trace.append({
            "n": n, "x_n": sign * x, "x_next": sign * x_new,
            "delta": delta, "residual": residual
        })
        if residual < eps:
            return sign * x_new, n, trace
        x = x_new
    raise RuntimeError(f"Точность не достигнута за {n_max} итераций.")


def ln_series(c: float, eps: float = EPS, n_max: int = N_MAX):
    """
    Вычисление ln(c) рядом:
    ln(c) = 2 * (y + y**3/3 + y**5/5 + ...), y=(c-1)/(c+1).
    Остановка: модуль очередного члена ряда меньше eps.
    """
    if c <= 0:
        raise ValueError("Аргумент натурального логарифма c должен быть положительным.")
    if eps <= 0:
        raise ValueError("Точность eps должна быть положительной.")
    y = (c - 1.0) / (c + 1.0)
    y2 = y * y
    power = y
    total = 0.0
    trace = []
    for k in range(n_max):
        term = 2.0 * power / (2 * k + 1)
        if abs(term) < eps:
            # Записываем последний проверенный член, который не прибавлялся к сумме.
            trace.append({
                "n": k + 1, "partial_sum": total,
                "term": term, "delta": abs(term), "added": False
            })
            return total, k, trace
        total += term
        trace.append({
            "n": k + 1, "partial_sum": total,
            "term": term, "delta": abs(term), "added": True
        })
        power *= y2
    raise RuntimeError(f"Точность не достигнута за {n_max} итераций.")


def print_root_trace(title, trace):
    print(f"\n{title}")
    print(f"{'n':>3} {'x_n':>18} {'|x_(n+1)-x_n|':>20} {'невязка':>16}")
    for row in trace:
        print(f"{row['n']:>3} {row['x_n']:>18.10f} "
              f"{row['delta']:>20.10e} {row['residual']:>16.10e}")


def print_ln_trace(trace):
    print("\nИтерации для ln(c)")
    print(f"{'n':>3} {'частичная сумма S_n':>24} {'|t_n|':>18} {'действие':>12}")
    for row in trace:
        action = "прибавлен" if row["added"] else "критерий стоп"
        print(f"{row['n']:>3} {row['partial_sum']:>24.12f} "
              f"{row['delta']:>18.10e} {action:>12}")


def main():
    sqrt_value, sqrt_n, sqrt_trace = sqrt_newton(A)
    cbrt_value, cbrt_n, cbrt_trace = cbrt_newton(B)
    ln_value, ln_n, ln_trace = ln_series(C)

    results = [
        ("sqrt(a)", A, sqrt_value, math.sqrt(A), sqrt_n),
        ("cbrt(b)", B, cbrt_value, B ** (1.0 / 3.0), cbrt_n),
        ("ln(c)", C, ln_value, math.log(C), ln_n),
    ]

    print("ЛР № 1 — вариант 19")
    print("Ситникова Ольга, группа БИ2-1")
    print(f"Параметры: a={A:g}, b={B:g}, c={C:g}, eps={EPS:g}")
    print(f"{'Функция':<12}{'Аргумент':>10}{'Результат':>18}"
          f"{'Эталон math':>18}{'Абс. ошибка':>16}{'Итераций':>10}")
    for name, arg, value, reference, n in results:
        print(f"{name:<12}{arg:>10.4g}{value:>18.10f}{reference:>18.10f}"
              f"{abs(value-reference):>16.6e}{n:>10}")

    print_root_trace("Итерации для sqrt(a)", sqrt_trace)
    print_root_trace("Итерации для cbrt(b)", cbrt_trace)
    print_ln_trace(ln_trace)


if __name__ == "__main__":
    main()
