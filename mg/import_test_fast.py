import imp
import json
import os
import pkgutil
import sys

database_file = os.path.join(os.path.dirname(__file__) or '.', 'import.cache')
module_db = {}

if os.path.exists(database_file):
  print "Database exists"
  with open(database_file, 'r') as fh:
    module_db = json.load(fh)

class CachingMetaImportFinder(object):
  def find_module(self, fullname, path=None):
    if '.' in fullname:
      # Only deal with top level packages
      return None
    if fullname in module_db:
      solution = module_db[fullname]
      if not solution:
        # Previously seen module that can't be cached
        print "Uncacheable", fullname
        return None
      print "Cached solution for", fullname, ":", solution
      fh = open(solution[1], solution[2][1]) if solution[0] else None
      return pkgutil.ImpLoader(fullname, fh, solution[1], solution[2])
    print "New module encountered:", fullname
    try:
      solution = imp.find_module(fullname)
    except:
      module_db[fullname] = None
      print "Uncacheable"
      return None
    module_db[fullname] = (solution[0] is not None, ) + solution[1:]
    return pkgutil.ImpLoader(fullname, *solution)

sys.meta_path.append(CachingMetaImportFinder())

import dials.util.version
print dials.util.version.dials_version()

with open(database_file, 'w') as fh:
  json.dump(module_db, fh, indent=2)
