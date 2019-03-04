from __future__ import print_function
import sys


def scrape_parameters(refinement_debug_log):
    parameters = []
    for record in open(refinement_debug_log):
        if "name mapping" in record:
            continue
        if "Parameter" in record and " : " in record:
            parameters.append(record.split()[-1])

    return parameters


def scrape_Detector(refinement_debug_log, parameters):
    values = {}
    sigmas = {}

    parameters = [p.replace("Detector1", "") for p in parameters]

    parameter = None
    for record in open(refinement_debug_log):
        if not record.strip():
            parameter = None
            continue
        if "Parameter Report" in record:
            continue
        if "Parameter " in record and record.strip().endswith(":"):
            name = record.split()[-1].replace(":", "")
            if name in parameters:
                parameter = name
        if "Value:" in record and parameter:
            values[parameter] = float(record.split()[-1])
        if "Sigma:" in record and parameter:
            sigmas[parameter] = float(record.split()[-1])

    return values, sigmas


refinement_debug_logs = sys.argv[1:]

parameters = scrape_parameters(refinement_debug_logs[0])

detector_parameters = []

for parameter in parameters:
    if parameter.startswith("Detector1"):
        detector_parameters.append(parameter.replace("Detector1", ""))

results = {}
for refinement_debug_log in refinement_debug_logs:
    values, sigmas = scrape_Detector(refinement_debug_log, parameters)
    results[refinement_debug_log] = values, sigmas

all_values = {}
for k in detector_parameters:
    all_values[k] = []
    for j, log in enumerate(refinement_debug_logs):
        all_values[k].append(results[log][0][k])

# convert milliradians to degrees
import math

mr2d = 180 / (1000 * math.pi)

for k in detector_parameters:
    values = all_values[k]
    mean = sum(values) / len(values)
    if "Tau" in k:
        diffs = [mr2d * (v - mean) for v in values]
    else:
        diffs = [(v - mean) for v in values]
    print("%14s" % k, end=" ")
    print(" ".join(["%6.3f" % d for d in diffs]))
