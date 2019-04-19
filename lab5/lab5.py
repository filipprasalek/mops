import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize_scalar
from sklearn.metrics import mean_squared_error


def parse_data(data):
    return [float(record.split()[1]) for record in data.read().split("\n") if record]


def read_input(raninfall_file, measurement_file):
    rainfall = open(raninfall_file, "r")
    measurement = open(measurement_file, "r")
    return parse_data(rainfall), parse_data(measurement)


def convolution_integral(tt, input, lam, t):
    result = 0
    for i in range(t):
        result += input[i] * 1 / tt * math.exp(-(t - i) / tt) * math.exp(-lam * (t - i))
    return result


def estimate_results(tt, rainfall, lam, t):
    estimated_result = []
    for i in range(t):
        estimated_result.append(convolution_integral(tt, rainfall, lam, i) if i > 160 else 0)
    return estimated_result


def calculate_mean_squared_error(tt, *args):
    input, output, lam, t = args
    return mean_squared_error(output, estimate_results(tt, input, lam, t))


lam = 0.004696
tt = 10
(rainfall, measurement) = read_input("opady", "dunaj")
t = len(rainfall)

estimated_result = estimate_results(tt, rainfall, lam, t)

actual_plot = plt.plot(range(t), measurement, label="Actual measurement")
estimated_plot = plt.plot(range(t), estimated_result, label="Estimated values")
plt.title('[Figure 1] Comparision between actual and estimated tracer concentration')
plt.xlabel('t [months]')
plt.ylabel('tracer concentration [TU]')
plt.legend()
plt.show()

best_tt = minimize_scalar(calculate_mean_squared_error, method='bounded', bounds=[1, 1000],
                          args=(rainfall, measurement, lam, t))
print(best_tt)