import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.timeseries import LombScargle
import os
import sys
from datetime import datetime
from dateutil import tz
import statistics
from scipy import signal
from random import random
import math

def moving_average_filter(t, y, m):
    """
    Arguments:
    t - time array
    y - value array
    m - kernel length

    Returns filtered value array
    """

    new_y = []

    for i, value in enumerate(y):
        sum = 0
        time_at_value = t[i]
        n = 0
        for j, time in enumerate(t):
            if (abs(time - time_at_value) <= m / 2):
                sum += y[j];
                n += 1
        sum /= n
        new_y.append(sum)

    return new_y

def low_pass_filter(t, y, f_max):
    return moving_average_filter(t, y, 1/f_max)

def high_pass_filter(t, y, f_min):
    low_pass = moving_average_filter(t, y, 1/f_min)
    return [y[i] - low_pass[i] for i in range(len(y))]

def band_pass_filter(t, y, f_min, f_max):
    return low_pass_filter(t, high_pass_filter(t, y, f_min), f_max)

def find_best_period(t, y, FAP=False, return_model=False):
    baseline = max(t) - min(t)
    ls = LombScargle(t, y)
    frequency, power = ls.autopower(minimum_frequency=1/baseline, maximum_frequency=1/2)
    periods = 1 / frequency
    best_power = power.max()
    best_period = periods[list(power).index(best_power)]
    if (FAP):
        FAP = ls.false_alarm_probability(best_power)
        if (return_model):
            return (best_period, FAP, ls)
        else:
            return (best_period, FAP)
    else:
        if (return_model):
            return (best_period, ls)
        else:
            return best_period

def main():

    # generate sparse, uneven times from 0-4000 days
    t = [i for i in range(4000) if random() > 0.80]

    sig = []
    for time in t:
        value = 0

        # add long-term magnetic cycles (3000 day period)
        value += math.sin(1/3000 * time * 2 * math.pi)

        # add rotational period (41 day period)
        value += 0.3 * math.sin(1/41 * time * 2 * math.pi)

        # add noise
        value += 3 * random()

        sig.append(value)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    ax1.plot(t, sig)
    ax1.set_title('Simulated S-values')
    ax1.axis([0, 4000, -2, 4])

    # filter out magnetic cycles (remove periods > 200 days)
    high_filtered = high_pass_filter(t, sig, 1/200)

    # filter out low periods (remove periods < 4 days)
    low_filtered = low_pass_filter(t, high_filtered, 1/4)

    ax2.plot(t, high_filtered)
    ax2.set_title('After 200-day high-pass filter')
    ax2.axis([0, 4000, -2, 4])

    ax3.plot(t, low_filtered)
    ax3.set_title('After additional 4-day low-pass filter')
    ax3.axis([0, 4000, -2, 4])
    ax3.set_xlabel('Time (days)')

    def print_best_period(t, y):
        best_period, fap = find_best_period(t, y, FAP=True)
        print("Period: {0:.2f}, FAP: {1:.4f}".format(best_period, fap))

    print_best_period(t, sig)
    print_best_period(t, high_filtered)
    print_best_period(t, low_filtered)

    plt.tight_layout()
    plt.show()

if (__name__ == '__main__'):
    main()
