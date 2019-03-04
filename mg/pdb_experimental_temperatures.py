from __future__ import division, absolute_import
from __future__ import print_function
from dials.util import debug_console
import Queue
import iotbx.cif
import gzip
import os
import random
import sys
import time
import threading
import traceback

# Disable output buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

pdb = "/dls/science/groups/scisoft/PDB-mmcif"
shadow = "/dls/tmp/wra62962/PDB"
try:
    os.mkdir(shadow)
except OSError:
    pass

start = time.time()
next_update = start + 15
incr_update = start + 1


def parse_temperature(string):
    if string in ("", "?", "None"):
        return None
    if string == "False":
        return False
    if "e" in string or "E" in string:
        sys.stderr.write("PARSE PROBLEM: " + string + "\n")
    return float(string)


def get_cached_temperature(filename):
    shadowfile = filename.replace(pdb, shadow)
    if os.path.exists(shadowfile):
        with open(shadowfile, "r") as fh:
            cached = fh.read()
        if cached:
            return parse_temperature(cached)
    temp = get_temperature(filename)
    try:
        os.mkdir(os.path.dirname(shadowfile))
    except OSError:
        pass
    with open(shadowfile, "w") as fh:
        fh.write(str(temp))
    return temp


def get_temperature(filename):
    with gzip.open(filename, "r") as fh:
        cif = iotbx.cif.reader(file_object=fh)
    model = cif.model()
    assert len(model) == 1, list(model)
    entry = model[list(model)[0]]
    if entry.get("_exptl.method") != "X-RAY DIFFRACTION":
        return False

    temp = entry.get("_diffrn.ambient_temp")
    if not temp:
        return None

    if isinstance(temp, basestring):
        return parse_temperature(temp)

    # scitbx array
    templist = filter(None, map(parse_temperature, temp))
    if not templist:
        return None
    return templist[0]


q_directories = Queue.Queue()
q_files = Queue.PriorityQueue()
q_results = Queue.Queue()

filecounter_lock = threading.Lock()
filecounter = 0

result_lock = threading.Lock()
resultcounter = 0
temperatures = []


def print_temperature_summary():
    total = len(temperatures)
    print("-" * 31)
    line = "{range:17s}:{count:6d} ({percentage:4.1f}%)"
    filtered = filter(None, temperatures)
    tempgiven = len(filtered)
    noxray = len(filter(lambda t: t is False, temperatures))
    notemp = total - tempgiven - noxray

    def write_line(rangestr, filterfn):
        c = len(filter(filterfn, filtered))
        if tempgiven > 0:
            percentage = 100 * c / tempgiven
        else:
            percentage = 0
        print(line.format(range=rangestr, count=c, percentage=percentage))

    write_line("       T <=   0K", lambda t: t <= 0 and t is not False)
    write_line("  0K < T <=  70K", lambda t: 0 < t <= 70)
    write_line(" 70K < T <=  95K", lambda t: 70 < t <= 95)
    write_line(" 95K < T <= 105K", lambda t: 95 < t <= 105)
    write_line("105K < T <= 160K", lambda t: 105 < t <= 160)
    write_line("160K < T <= 265K", lambda t: 160 < t <= 265)
    write_line("265K < T <= 295K", lambda t: 265 < t <= 295)
    write_line("295K < T        ", lambda t: 295 < t)
    print("{range:^17s}:{count:6d}".format(range="T not given", count=notemp))
    print("{range:^17s}:{count:6d}".format(range="not xray", count=noxray))
    print("{range:^17s}:{count:6d}".format(range="total", count=total))
    print("-" * 31)


class Worker(threading.Thread):
    """ Thread doing things in the background """

    def __init__(self):
        threading.Thread.__init__(self)
        self.daemon = True
        self.start()


class DirectoryIndexer(Worker):
    """ Thread finding files in directories from q_directories into q_files """

    def run(self):
        global filecounter
        global filecounter_lock
        while True:
            dirname = q_directories.get()
            try:
                files = filter(
                    os.path.isfile,
                    (os.path.join(dirname, f) for f in os.listdir(dirname)),
                )
                with filecounter_lock:
                    filecounter = filecounter + len(files)
                for f in files:
                    q_files.put((random.randrange(1000000), f))
            except Exception:
                # An exception happened in this thread
                traceback.print_exc()
            finally:
                # Mark this task as done, whether an exception happened or not
                q_directories.task_done()


class FileParser(Worker):
    """ Thread parsing from q_files into q_results """

    def run(self):
        while True:
            _, filename = q_files.get()
            try:
                q_results.put(get_cached_temperature(filename))
            except Exception:
                # An exception happened in this thread
                print("\n\nFailed to read temperature from", filename)
                traceback.print_exc()
            finally:
                # Mark this task as done, whether an exception happened or not
                q_files.task_done()


class ResultParser(Worker):
    """ Thread aggregating from q_results """

    def run(self):
        global result_lock
        global resultcounter
        global temperatures
        global next_update
        global incr_update

        while True:
            temperature = q_results.get()
            try:
                with result_lock:
                    resultcounter = resultcounter + 1
                    temperatures.append(temperature)
                    if time.time() > next_update:
                        print_temperature_summary()
                        next_update = time.time() + 15

                if time.time() > incr_update:
                    now = time.time()
                    incr_update = now + 1
                    completed = resultcounter / filecounter
                    percentage = completed * 100
                    eta = start + ((now - start) / completed) - now
                    if eta > 120:
                        timestr = "{eta:1.0f} minutes".format(eta=eta / 60)
                    else:
                        timestr = "{eta:1.0f} seconds".format(eta=eta)
                    print(
                        "{percentage:4.1f}% - ETA: {timestring} - {done} of {total} files parsed".format(
                            done=resultcounter,
                            total=filecounter,
                            timestring=timestr,
                            percentage=percentage,
                        )
                    )

            except Exception:
                # An exception happened in this thread
                traceback.print_exc()
            finally:
                # Mark this task as done, whether an exception happened or not
                q_results.task_done()


try:
    print("Starting DirectoryIndexer")
    for x in range(5):
        DirectoryIndexer()
    print("Loading base directory")
    directories = filter(os.path.isdir, (os.path.join(pdb, f) for f in os.listdir(pdb)))
    for d in directories:
        q_directories.put(d)
    # q_directories.put(directories[180])
    print("Starting CIF parsers")
    for x in range(3):
        FileParser()
    print("Starting result parser")
    ResultParser()
    print("Running normally")

    q_directories.join()
    print("All directories read")
    q_files.join()
    print("All files parsed")
    q_results.join()
    print("All results aggregated")
    print_temperature_summary()
except KeyboardInterrupt:
    print("Killed by Ctrl+C")
    sys.exit(1)
