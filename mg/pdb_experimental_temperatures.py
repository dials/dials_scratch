from __future__ import absolute_import, division, print_function

import collections
import concurrent.futures
import gzip
import iotbx.cif
import itertools
import json
import os
import sys
import time
import threading
from tqdm import tqdm

# Disable output buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), "w", 0)

pdb = "/dls/science/groups/scisoft/PDB-mmcif"
dbjson = "/dls/tmp/wra62962/PDB.json"
try:
    with open(dbjson, "r") as fh:
        db = json.load(fh)
except IOError:
    db = {}
db_lock = threading.Lock()


def parse_temperature(string):
    if string in ("", "?", "None"):
        return None
    if string == "False":
        return False
    if "e" in string or "E" in string:
        sys.stderr.write("PARSE PROBLEM: " + string + "\n")
    return float(string)


def get_year_temperature(filename):
    with gzip.open(filename, "r") as fh:
        cif = iotbx.cif.reader(file_object=fh)
    model = cif.model()
    assert len(model) == 1, list(model)
    entry = model[list(model)[0]]
    if entry.get("_exptl.method") != "X-RAY DIFFRACTION":
        return False, False
    year = entry.get("_pdbx_database_status.recvd_initial_deposition_date")
    if year and "-" in year:
        year = int(year.split("-")[0])
    else:
        year = False

    temp = entry.get("_diffrn.ambient_temp")
    if not temp:
        return year, None

    if isinstance(temp, basestring):
        return year, parse_temperature(temp)

    # scitbx array
    templist = filter(None, map(parse_temperature, temp))
    if not templist:
        return year, None
    return year, templist[0]


def print_temperature_summary(temperature_dict):
    for year in sorted(temperature_dict):
        temperatures = temperature_dict[year]
        total = len(temperatures)
        print("-" * 31)
        print("{year:^17s}:{yearno:6d}".format(year="Year", yearno=year))
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


def read_directory(dirname):
    """Function finding files in directories"""

    files = filter(
        os.path.isfile, (os.path.join(dirname, f) for f in os.listdir(dirname))
    )
    return tuple(files)


def read_file(filename):
    pdb_code = os.path.basename(filename)
    if pdb_code not in db:
        year_temp = get_year_temperature(filename)
        with db_lock:
            db[pdb_code] = year_temp
    return db[pdb_code]


def find_all_files():
    print("Reading directories directory...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as directory_reader:
        print("Reading top-level directory...")
        directories = list(os.listdir(pdb))
        all_directories = []
        with tqdm(total=len(directories)) as pbar:
            for d in directories:
                full_directory_path = os.path.join(pdb, d)
                pbar.update()
                if not os.path.isdir(full_directory_path):
                    continue
                all_directories.append(
                    directory_reader.submit(read_directory, full_directory_path)
                )
        print("Done. Reading subdirectories...")
        with tqdm(total=len(all_directories)) as pbar:
            for result in concurrent.futures.as_completed(all_directories):
                pbar.update()
    print("Done.")
    return list(itertools.chain.from_iterable(d.result() for d in all_directories))


files = find_all_files()
num_files = len(files)

with concurrent.futures.ThreadPoolExecutor(max_workers=10) as file_reader:
    resultcounter = 0
    temperatures = collections.defaultdict(list)
    start = time.time()
    next_update = start + 15
    active_futures = set()
    print("Creating 50 futures")
    with tqdm(total=50) as pbar:
        for n in range(50):
            if files:
                active_futures.add(file_reader.submit(read_file, files.pop()))
            pbar.update()

    print("Reading %d files" % num_files)
    try:
        with tqdm(total=num_files) as pbar:
            while active_futures:
                for result in concurrent.futures.as_completed(active_futures):
                    pbar.update()
                    year, temperature = result.result()
                    resultcounter = resultcounter + 1
                    active_futures.remove(result)
                    if files:
                        active_futures.add(file_reader.submit(read_file, files.pop()))
                    if year is False and temperature is False:
                        continue

                    temperatures[year].append(temperature)
                    if time.time() > next_update:
                        print("\n")
                        print_temperature_summary(temperatures)
                        print("")
                        with db_lock:
                            json_string = json.dumps(db)
                        with open(dbjson, "w") as fh:
                            fh.write(json_string)
                        json_string = None
                        next_update = time.time() + 60

        print("All results aggregated")
        print_temperature_summary(temperatures)
        with db_lock:
            with open(dbjson, "w") as fh:
                json.dump(db, fh)
    except KeyboardInterrupt:
        print("Cancelling processes...")
        for f in active_futures:
            f.cancel()
        print("Saving progress...")
        with db_lock:
            with open(dbjson, "w") as fh:
                json.dump(db, fh)
        print("Waiting for processes to finish...")
        sys.exit(1)
