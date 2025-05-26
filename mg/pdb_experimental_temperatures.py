import collections
import concurrent.futures
from gemmi import cif
import gzip
import itertools
import json
import os
import sys
import time
from tqdm import tqdm

pdb = "/dls/science/groups/scisoft/PDB-mmcif"
dbjson = "/dls/tmp/wra62962/PDB.json"
if not os.path.exists(pdb):
    pdb = "/home/markus/pdb/pdb"
    dbjson = "/home/markus/pdb/state.json"
csv = dbjson[:-4] + "csv"

try:
    with open(dbjson, "r") as fh:
        db = json.load(fh)
except IOError:
    db = {}


def parse_temperature(string):
    if string in ("", "?", "None", "."):
        return None
    if string == "False":
        return False
    if "e" in string or "E" in string:
        sys.stderr.write("PARSE PROBLEM: " + string + "\n")
    try:
        return float(string)
    except ValueError:
        sys.stderr.write("Can not parse temperature: " + string + "\n")
        raise
        return False


def cif_get(block, item):
    column = block.find_values(item)
    if len(column) == 0:
        return False
    if len(column) == 1:
        return column.str(0)
    for value in column:
        if value:
            return value
    return False


def get_year_temperature(filename):
    doc = cif.read(filename)
    block = doc.sole_block()

    year = cif_get(block, "_pdbx_database_status.recvd_initial_deposition_date")
    if year:
        if "-" in year:
            year = int(year.split("-")[0])
        else:
            year = False
    if cif_get(block, "_exptl.method") != "X-RAY DIFFRACTION":
        return year, False
    temp = cif_get(block, "_diffrn.ambient_temp")
    if not temp:
        return year, None

    return year, parse_temperature(temp)


def print_temperature_summary(temperature_dict):
    all_years = sorted(temperature_dict)
    years = all_years[0:8]
    print("-" * (17 + 15 * len(years)))
    while years:
        all_years = all_years[8:]
        print(
            "{:^17s}:".format("Year")
            + "         ".join("{:6d}".format(year) for year in years)
        )

        total = {year: len(temperature_dict[year]) for year in years}
        noxray = {
            year: len([t for t in temperature_dict[year] if t is False])
            for year in years
        }
        filtered = {
            year: [t for t in temperature_dict[year] if t not in (False, None)]
            for year in years
        }
        tempgiven = {year: len(filtered[year]) for year in years}
        notemp = {year: total[year] - tempgiven[year] - noxray[year] for year in years}

        line = "{:17s}:" + " ".join("{:6d} ({:4.1f}%)" for year in years)

        def write_line(rangestr, filterfn):
            args = [rangestr]
            for year in years:
                c = len([item for item in filtered[year] if filterfn(item)])
                if tempgiven[year] > 0:
                    percentage = 100 * c / tempgiven[year]
                else:
                    percentage = 0
                if percentage > 99.9:
                    percentage = 99.9
                args.append(c)
                args.append(percentage)
            print(line.format(*args))

        write_line("       T <=   0K", lambda t: t <= 0 and t is not False)
        write_line("  0K < T <=  70K", lambda t: 0 < t <= 70)
        write_line(" 70K < T <=  95K", lambda t: 70 < t <= 95)
        write_line(" 95K < T <= 105K", lambda t: 95 < t <= 105)
        write_line("105K < T <= 160K", lambda t: 105 < t <= 160)
        write_line("160K < T <= 265K", lambda t: 160 < t <= 265)
        write_line("265K < T <= 295K", lambda t: 265 < t <= 295)
        write_line("295K < T        ", lambda t: 295 < t)

        line = "{:17s}:" + "         ".join("{:6d}" for year in years)
        args = ["T not given"] + [notemp[y] for y in years]
        print(line.format(*args))
        args = ["not xray"] + [noxray[y] for y in years]
        print(line.format(*args))
        args = ["total"] + [total[y] for y in years]
        print(line.format(*args))
        print("-" * (17 + 15 * len(years)))
        years = all_years[0:8]
    write_csv(temperature_dict)


def write_csv(temperature_dict):
    with open(csv, "w") as fh:
        all_years = sorted(temperature_dict)
        fh.write("Category," + ",".join([str(year) for year in all_years]) + "\n")

        filtered = {
            year: [t for t in temperature_dict[year] if t not in (False, None)]
            for year in all_years
        }

        def write_line(rangestr, filterfn):
            row = ["'" + rangestr + "'"]
            for year in all_years:
                c = len([item for item in filtered[year] if filterfn(item)])
                row.append(c)
            fh.write(",".join(str(c) for c in row) + "\n")

        write_line("       T <=   0K", lambda t: t <= 0 and t is not False)
        write_line("  0K < T <=  70K", lambda t: 0 < t <= 70)
        write_line(" 70K < T <=  95K", lambda t: 70 < t <= 95)
        write_line(" 95K < T <= 105K", lambda t: 95 < t <= 105)
        write_line("105K < T <= 160K", lambda t: 105 < t <= 160)
        write_line("160K < T <= 265K", lambda t: 160 < t <= 265)
        write_line("265K < T <= 295K", lambda t: 265 < t <= 295)
        write_line("295K < T        ", lambda t: 295 < t)

        total = {year: len(temperature_dict[year]) for year in all_years}
        noxray = {
            year: len([t for t in temperature_dict[year] if t is False])
            for year in all_years
        }
        tempgiven = {year: len(filtered[year]) for year in all_years}
        notemp = {
            year: total[year] - tempgiven[year] - noxray[year] for year in all_years
        }
        fh.write("'T not given'," + ",".join(str(notemp[y]) for y in all_years) + "\n")
        fh.write("'not xray'," + ",".join(str(noxray[y]) for y in all_years) + "\n")
        fh.write("total," + ",".join(str(total[y]) for y in all_years) + "\n")


def read_directory(dirname):
    """Function finding files in directories"""

    files = filter(
        os.path.isfile, (os.path.join(dirname, f) for f in os.listdir(dirname))
    )
    return tuple(files)


def read_file(filename):
    pdb_code = os.path.basename(filename)
    if pdb_code in db:
        return pdb_code, db[pdb_code]
    else:
        return pdb_code, get_year_temperature(filename)


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

with concurrent.futures.ProcessPoolExecutor(max_workers=8) as file_reader:
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
                    pdb_code, y_t = result.result()
                    db[pdb_code] = y_t
                    year, temperature = y_t
                    resultcounter = resultcounter + 1
                    active_futures.remove(result)
                    if files:
                        active_futures.add(file_reader.submit(read_file, files.pop()))
                    if year is False and temperature is False:
                        continue

                    temperatures[year].append(temperature)
                    if time.time() > next_update:
                        print("\n")
                        json_string = json.dumps(db)
                        with open(dbjson, "w") as fh:
                            fh.write(json_string)
                        json_string = None
                        print_temperature_summary(temperatures)
                        print("")
                        next_update = time.time() + 60

        with open(dbjson, "w") as fh:
            json.dump(db, fh)
        print("All results aggregated")
        print_temperature_summary(temperatures)
    except KeyboardInterrupt:
        print("Cancelling processes...")
        for f in active_futures:
            f.cancel()
        print("Saving progress...")
        with open(dbjson, "w") as fh:
            json.dump(db, fh)
        print("Waiting for processes to finish...")
        sys.exit(1)
