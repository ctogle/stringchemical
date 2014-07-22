
import stringchemical as sc

if hasattr(sc, 'simulate'):
	print 'stringchemical module imported'

import chemical.scripts.chemicallite as lcl

if hasattr(lcl, 'module_name'):
	print 'chemicallite module imported as', lcl.module_name
