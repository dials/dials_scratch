. ~/svn/setup.sh

export P6M_60_PANEL=1

mkdir ${1}
cd ${1}
dials.index ../datablock-${1}.json ../strong-${1}.pickle space_group=P4 \
reflections_per_degree=None
dials.refine indexed.pickle experiments.json hierarchy_level=1 \
reflections_per_degree=None ../refine.phil
