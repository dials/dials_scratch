. ~/svn/setup.sh

run=$1

export P6M_60_PANEL=1

dials.import /dls/mx-scratch/gw56/metrology_paper/i03/data0/${run}
dials.find_spots datablock.json nproc=20 shoeboxes=false
dials.index datablock.json strong.pickle space_group=P4 \
reflections_per_degree=None
dials.refine indexed.pickle experiments.json hierarchy_level=1 \
reflections_per_degree=None ../refine.phil