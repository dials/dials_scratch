. ~/svn/setup.sh

export P6M_60_PANEL=1

run=${1}
data=/dls/i03/data/2016/cm14451-4/processing/gw/metro/1/data/${run}

dials.import template=${data}/INS3_1_####.cbf
dials.find_spots datablock.json nproc=20 shoebox=false
dials.index datablock.json strong.pickle space_group=I213 \
reflections_per_degree=None

