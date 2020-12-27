from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import posixpath
import sys
import time
import signal
from math import *
from scipy.signal import find_peaks
import wfdb
import argparse

parser = argparse.ArgumentParser(description="Takes a record name, which is in RECORDS.txt, and parses it to detect apnea by given timestamp")
parser.add_argument("record", help="name of the record to use", type=str)
parser.add_argument("timestamp", help="timestamp to use", type=int)
parser.add_argument("-p", "--plot", help="display each timestamp of record", default=0,
                    action="count")
args = parser.parse_args()

samp = 0
coord = np.array([])
timestamp = int(args.timestamp)

tmp_record = wfdb.rdrecord(args.record, channels=[0], pn_dir='apnea-ecg')
step = int(tmp_record.fs)
previousRMSSD = 0
apnea = ""
tour=0

if timestamp == 0:
    timestamp = tmp_record.sig_len - 1

try:
    while (samp < tmp_record.sig_len):
        sampto = step * timestamp + samp
        if (step * timestamp + samp >= tmp_record.sig_len):
            sampto = tmp_record.sig_len - 1
        record = wfdb.rdrecord(args.record, pn_dir='apnea-ecg', channels=[0], sampfrom = samp, sampto = sampto)
        ann = wfdb.rdann(args.record, 'qrs', sampfrom = samp, sampto = step * timestamp + samp, shift_samps=True, pn_dir='apnea-ecg')
        intervalRR = np.delete(ann.sample, 0)
        for i in range(1, ann.sample.size):
            intervalRR[i - 1] = (ann.sample[i] - ann.sample[i-1])
        RMSSD = 0
        for i in range(0, intervalRR.size):
            RMSSD += intervalRR[i]*intervalRR[i]
        RMSSD = sqrt(RMSSD / ann.sample.size - 1)
        if (previousRMSSD != 0):
            if (RMSSD > 89 and previousRMSSD > 87):
                apnea = "A"
            else:
                apnea = "N"
        print(time.strftime('%H:%M:%S', time.gmtime(tour * timestamp)) + ": ", RMSSD, apnea)
        x = np.delete(ann.sample, 0)
        for i in range(0, x.size):
            x[i] = x[i] / 100
        if args.plot >= 1:
            plt.plot(x, intervalRR)
            plt.title('Record of Interval from ' + args.record + ' (' + time.strftime('%H:%M:%S', time.gmtime(tour * timestamp)) + ')')
            plt.xlabel('Time (s)')
            plt.ylabel('Interval (ms)')
            plt.show()
        samp += timestamp * step
        tour+=1
        previousRMSSD = RMSSD
except KeyboardInterrupt:
    print()