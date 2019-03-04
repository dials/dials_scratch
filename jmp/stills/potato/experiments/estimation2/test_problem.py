from dials_scratch.jmp.stills.potato.profile_refiner import ProfileRefiner
import pickle
import sys

data = pickle.load(open(sys.argv[1]))
s0 = data["s0"]
s2_list = data["s2_list"]
ctot_list = data["ctot_list"]
Sobs_list = data["Sobs_list"]


refiner = ProfileRefiner(s0, s2_list, ctot_list, Sobs_list)
refiner.refine()
params = refiner.parameters
