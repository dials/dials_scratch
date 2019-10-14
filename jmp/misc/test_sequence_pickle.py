from dials.model.serialize import load
import pickle

sequence = load.sequence(
    "/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/sequence.json"
)

text = pickle.dumps(sequence)

sequence2 = pickle.loads(text)

sequence2.get_detector()
