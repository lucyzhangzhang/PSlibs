#!/usr/bin/python3

import urllib.request
import urllib.parse

# Arabidopsis = 3702
#	Core organisms, Blasted agains the NCBI database
# Eutrema = 72664
#	Periphery organism, only Blasted against core organisms

# test
url = 'http://string-db.org/api/tsv/network?identifier=IPS1&species=3702'
f = urllib.request.urlopen(url)
print(f.read().decode('utf-8'))
