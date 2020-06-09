# Object Store Prototyping

This document describes work on the object store implementation and prototyping. It summarizes the existing system and describes the differences with the proposed new structure. It then goes on to describe how to use the prototype tools to run the system yourself.

## Existing Infrastructure

The existing infrastructure is based on a service called TOM, which writes each image to a separate file in GPFS. This is only used for grid scans, not the automated analysis pipeline - because the latter can run on complete data sets after collection, but this may be a long time before you start getting feedback.

The workflow is triggered by GDA. Once started, for every image read from the detector there is several services that interact:

![TOM Sequence](http://www.plantuml.com/plantuml/proxy?src=https://raw.githubusercontent.com/dials/dials_scratch/master/gw/object_store/doc/tom.puml)

- The **DAQ** sits on the beamline. The DAQ tells the detector to start collecting, waits for messages from the detector containing image data, collects them and directs the data to the FileWriter. The DAQ _also_ forwards the detector image messages to the TOM service
- **TOM** listens for image messages from the DAQ. TOM writes each of these as a separate file to a temporary directory in the visit folder on **GPFS**.
- When the workflow is triggered, the **FileWatcher** is set up with the known
  filesystem location that TOM writes it files to. It checks for each of these
  in turn, intermittently polling the filesystem for the existence of the
  file. When it detects the next image is present on disk, it sends a message to the spotfinding service telling it the filename and index of the new image
- The **Spotfinder** service, once notified, grabs the image from disk, runs the analysis and sends the results onwards, back to GDA

### Drawbacks
There are a few issues with using the filesystem as the intermediatary here, which are mainly consequences of abstracting the storage to a POSIX-compatible filesystem:
- Large number of files being pushed onto filesystem storage. This causes high load on the service and enumerating the contents of large directories is very slow
- Response time: GPFS can sometimes be unpredictable in when filesystem changes become available to different clients. e.g. TOM has written a file and the FileWatcher can see it exists, but the Spotfinder can't yet. There are deliberate delays built into the process to avoid this, and if the file doesn't exist services wait and then try again some short period later. This causes overall lag in the end-to-end response time

## Proposed Infrastructure

![TOM Sequence](http://www.plantuml.com/plantuml/proxy?src=https://raw.githubusercontent.com/dials/dials_scratch/master/gw/object_store/doc/proposed.puml)

- The **DAQ** behaves the same as before; communicates with the detector, sends the data to the FileWriter (and TOM while prototyping)
- A new **CephWriter** service waits for images from the DAQ, pushes them into the **CEPH** object store and then sends a message to **Zocalo** notifying it that the image is available
- The first available **Spotfinder** gets _immediately_ notified, fetches the image from the object store and then does the computation. The results is then returned to GDA, or whatever observing service triggered the workflow.

## Benefits
- Data written to the object store is immediately available to other readers once it is completed
- Objects in the object store can be set to expire, so no manual management of data removal is required
- There is no polling or waiting on filesystems

## CEPH Prototype
Scientific Computing has kindly supplied a test server running CEPH. It runs on RamDisk, so is constrained in the amount of data it can store. However, this is good as a test. The trial server is currently running on  `cs04r-sc-com99-25.diamond.ac.uk:7480`

## Running the trial
- Have `dials_scratch` and `dlstbx` installed into your dials environment
- Make sure you have the minio api - `libtbx.pip install minio`
- Save the Secret Key to an environment variable - for the internal testing ceph credentials, see the "Re: RAM-based ceph Instance" email or ask directly.
```
export SECRET_KEY=<value>
```

### Loading data (skip if using existing)
Load the data you want to analyse into ceph with the form:
```
dials.python dials_scratch/gw/object_store/ingest.py --host <HOST> <FILE_master.h5> <DCID>
```
Where `<HOST>` is a URL to the S3 endpoint, including access credentials e.g. 
`s3://access_key:secret_key@hostname:port`.


For testing I have pre-loaded the "large test" Neil Patterson high-res Proteinase-K dataset from https://ispyb.diamond.ac.uk/dc/visit/cm23003-2/id/3474639 with:
```
dials.python dials_scratch/gw/object_store/ingest.py
 --host "s3://XAVGY5RBY28BTXRENEII:${SECRET_KEY}@cs04r-sc-com99-25.diamond.ac.uk:7480" \
  /dls/i03/data/2019/cm23003-2/proteinasek/protk_1/protk_1_1_master.h5 \
  3474639
 ```

 ### Create Zocalo server and services
 Create a test activeMQ server to run zocalo on; These examples use `ZOCHOST` as a placeholder for whichever machine you are running activeMQ on:
 ```
 module load activemq/zocdev
```
Now _on the same machine_ run the status and queue monitors and launch the dispatcher and CephWatcher services - for convenience, in a single terminal with tmux:
```
tmux new-session \; \
  send-keys "clear && dlstbx.status_monitor --test" C-m \; \
  split-window -v \; \
  send-keys "sleep 5 && clear && dlstbx.queue_monitor --test" C-m \; \
  split-window -h \; \
  send-keys "clear && dlstbx.service -s DLSDispatcher --test" C-m \; \
  split-window -v \; \
  send-keys "clear && dlstbx.service -s DLSCephWatcher --test" C-m \;
```
(if launching on a separate machine you will need to pass the `--stomp-host=ZOCHOST` argument). You can kill these processes with `Ctrl-C` and `Ctrl-D`, which will close each window in turn. Or if you prefer not to use tmux, separately:
```
dlstbx.service -s DLSDispatcher --test
dlstbx.service -s DLSCephWatcher --test
dlstbx.status_monitor --test
dlstbx.queue_monitor --test
```

The **CephWatcher** service is to effiently inject a whole set of object-store-images to Per-Image-Analysis. Without this, it's hard to get a large flood of requests into Zocalo to test throughput.

### PIA Services

Create a `run_pia.sh` script to load your environment and start a PIA service:
```
#!/bin/bash
#$ -j y
. /path/to/dials/build/setpaths.sh
dlstbx.service -s DLSPerImageAnalysis --test --stomp-host=ZOCHOST
```
Now you can launch as many as you want with qsub - `1-20` launches 20 processes:
```
module load global/testcluster
qsub -N pia -t 1-20 -pe smp 4 ./run_pia.sh
```
(each PIA has been set to use 4 slots, to prevent too many running on the same machine for these tests)

### Submission

Use the submit script in `dials_scratch/gw/object_store/submit.py` to submit the endpoint (with bucket name) like this:
```
dials.python dials_scratch/gw/object_store/submit.py \
  "s3://XAVGY5RBY28BTXRENEII:${SECRET_KEY}@cs04r-sc-com99-25.diamond.ac.uk:7480/3474639'
```
This will send a message, then wait for the responses, printing them out as they come.