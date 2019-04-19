import matplotlib.pyplot as plt
import math


def parse_data(data):
    return [float(record.split()[1]) for record in data.read().split("\n") if record]


def read_input(raninfall_file, measurement_file):
    rainfall = open(raninfall_file, "r")
    measurement = open(measurement_file, "r")
    return parse_data(rainfall), parse_data(measurement)


def convolution_integral(input, lam, t, tt):
    result = 0
    for i in range(t):
        result += input[i] * 1 / tt * math.exp(-(t - i) / tt) * math.exp(-lam * (t - i))
    return result


lam = 0.004696
tt = 10
(rainfall, measurement) = read_input("opady", "dunaj")
t = len(rainfall)

estimated_result = []
for i in range(t):
    estimated_result.append(convolution_integral(rainfall, lam, i, tt) if i > 160 else 0)

actual_plot = plt.plot(range(t), measurement, label="Actual measurement")
estimated_plot = plt.plot(range(t), estimated_result, label="Estimated values")
plt.title('[Figure 1] Comparision between actual and estimated tracer concentration')
plt.xlabel('t [months]')
plt.ylabel('tracer concentration [TU]')
plt.legend()
plt.show()

