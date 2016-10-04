. ~/svn/setup.sh

export P6M_60_PANEL=1
suffix=$1
cd ${suffix}
dials.refine ../reflections_${suffix}.pickle ../experiments_${suffix}.json \
hierarchy_level=1 reflections_per_degree=None ../refine.phil

