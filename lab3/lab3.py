import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve


def calculate_temp_no_atmosphere(s, a, sigma):
    nominator = s * (1 - a)
    denominator = 4 * sigma
    return math.pow(nominator // denominator, 1 / 4)


def balance_equation(vars, *constants):
    T_a, T_s = vars
    t_a, a_a, a_s, t_a_prim, a_a_prim, c, s, sigma = constants
    eq_1 = -t_a * (1 - a_s) * s / 4 + c * (T_s - T_a) + sigma * math.pow(T_s, 4) * (1 - a_a_prim) - sigma * math.pow(
        T_a, 4)
    eq_2 = - (1 - a_a - t_a + a_s * t_a) * s / 4 - c * (T_s - T_a) - sigma * math.pow(T_s, 4) * (
            1 - t_a_prim - a_a_prim) + 2 * sigma * math.pow(T_a, 4)
    return eq_1, eq_2


def calculate_temp_with_atmosphere(t_a, a_a, a_s, t_a_prim, a_a_prim, c, s, sigma):
    return fsolve(balance_equation, np.array([0, 0]), args=(t_a, a_a, a_s, t_a_prim, a_a_prim, c, s, sigma))


solar_const = 1366
mean_albedo = 0.3
stefan_boltzmann_const = 5.67 * math.pow(10, -8)

temp_no_atmosphere = calculate_temp_no_atmosphere(s=solar_const, a=mean_albedo, sigma=stefan_boltzmann_const)
print("Earth mean temperature, assuming no atmosphere: %10.4f K" % temp_no_atmosphere)

short_wave_trans = 0.53
long_wave_trans = 0.06
surface_albedo_short = 0.19
atmosphere_albedo_short = 0.30
atmosphere_albedo_long = 0.31
c = 2.7

temp_with_atmosphere = calculate_temp_with_atmosphere(t_a=short_wave_trans,
                                                      a_a=atmosphere_albedo_short,
                                                      a_s=surface_albedo_short,
                                                      t_a_prim=long_wave_trans,
                                                      a_a_prim=atmosphere_albedo_long,
                                                      c=c,
                                                      s=solar_const,
                                                      sigma=stefan_boltzmann_const
                                                      )

print("Earth atmosphere mean temperature, considering atmosphere: %10.4f K" % temp_with_atmosphere[0])
print("Earth surface mean temperature, considering atmosphere: %10.4f K" % temp_with_atmosphere[1])

solar_const_range = range(math.floor(0.8 * solar_const), math.floor(1.2 * solar_const))
surface_temperatures = []
atmosphere_temperatures = []
for i in solar_const_range:
    temperature_pair = calculate_temp_with_atmosphere(t_a=short_wave_trans,
                                                      a_a=atmosphere_albedo_short,
                                                      a_s=surface_albedo_short,
                                                      t_a_prim=long_wave_trans,
                                                      a_a_prim=atmosphere_albedo_long,
                                                      c=c,
                                                      s=i,
                                                      sigma=stefan_boltzmann_const
                                                      )
    atmosphere_temperatures.append(temperature_pair[0])
    surface_temperatures.append(temperature_pair[1])

plt.plot(solar_const_range, atmosphere_temperatures)
plt.title('[Figure 1] Relationship between solar constant and mean atmosphere temperature')
plt.xlabel('Solar constant [W/m^2]')
plt.ylabel('Mean atmosphere temperature [K]')
plt.show()

plt.plot(solar_const_range, surface_temperatures)
plt.title('[Figure 2] Relationship between solar constant and mean surface temperature')
plt.xlabel('Solar constant [W/m^2]')
plt.ylabel('Mean surface temperature [K]')
plt.show()

albedo_glacial = surface_albedo_short
atmosphere_temperatures = []
surface_temperatures = []
for i in solar_const_range:
    temperature_pair = calculate_temp_with_atmosphere(t_a=short_wave_trans,
                                                      a_a=atmosphere_albedo_short,
                                                      a_s=albedo_glacial,
                                                      t_a_prim=long_wave_trans,
                                                      a_a_prim=atmosphere_albedo_long,
                                                      c=c,
                                                      s=i,
                                                      sigma=stefan_boltzmann_const
                                                      )
    atmosphere_temperatures.append(temperature_pair[0])
    surface_temperatures.append(temperature_pair[1])
    albedo_glacial = 0.62 if temperature_pair[1] < 273 else 0.19

plt.plot(solar_const_range, atmosphere_temperatures)
plt.title('[Figure 3] Relationship between solar constant and mean atmosphere temperature with glacial mechanism')
plt.xlabel('Solar constant [W/m^2]')
plt.ylabel('Mean atmosphere temperature [K]')
plt.show()

plt.plot(solar_const_range, surface_temperatures)
plt.title('[Figure 4] Relationship between solar constant and mean surface temperature with glacial mechanism')
plt.xlabel('Solar constant [W/m^2]')
plt.ylabel('Mean surface temperature [K]')
plt.show()