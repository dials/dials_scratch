
def make_plot(filename, title):

  with open(filename) as infile:
    from collections import defaultdict
    class Data(object):
      def __init__(self):
        self.x = []
        self.y = []
    data = defaultdict(Data)
    for line in infile.readlines():
      tokens = line.split()
      dset = int(tokens[0])
      num = int(tokens[1])
      sigma = float(tokens[2])
      if dset == 0:
        print line
      data[dset].x.append(num)
      data[dset].y.append(sigma)

    from matplotlib import pylab
    max_x = 2.0

    fig = pylab.figure(figsize=(6,6/1.6), dpi=300)
    for key in sorted(data.keys()):
      x = data[key].x
      y = data[key].y
      l = str(key)
      x = [2.0 * xx / len(x) for xx in x]
      pylab.plot(x, y, label=l)
    pylab.ylim((0,0.6))
    pylab.legend(loc="best")
    pylab.title(title)
    pylab.xlabel("Oscillation (degrees)")
    pylab.ylabel("Sigma_M")
    #pylab.show()

    from os.path import splitext
    outfile = "%s.png" % splitext(filename)[0]
    fig.savefig(outfile, dpi=300, bbox_inches="tight")


make_plot("results/mosaicity_old.txt", "ln(L) = sum(ln(P_i))")
make_plot("results/mosaicity_new.txt", "ln(L) = sum(n_i * ln(P_i)) - n_tot * ln(P_tot)")
make_plot("results/mosaicity_new2.txt", "ln(L) = sum(n_i * ln(P_i))")
make_plot("results/mosaicity_new3.txt", "ln(L) = sum(n_i * ln(P_i)) - n_tot * P_tot")
