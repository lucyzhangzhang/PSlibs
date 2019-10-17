#!/usr/bin/env python3
# need to have Intermine installed
# you can do a local install through conda

from intermine.webservice import Service
import argparse
import os

usage = "usage %prog GeneList \"OrganismShortName\""

parser = argparse.ArgumentParser(description = "Interacts with the Phytozome API")
parser.add_argument("-l", "--list", dest = "geneList", required = True, 
        metavar = "Gene List", help = "List of names of genes you want to search for, each on their own line")
parser.add_argument("-n", "--name", dest = "Sname", required = True, type = str,
        metavar = "Organism shortname (i.e. A. thaliana)", help = "Name of the organism that the gene comes from")

service = Service("https://phytozome.jgi.doe.gov/phytomine/service")

args = parser.parse_args()

#path of the input file
fpath = os.path.dirname(args.geneList)

with open(args.geneList) as gList:
    for gName in gList:
        gName = gName.rstrip()
        query = service.new_query("Homolog")
        query.add_view(
            "gene.primaryIdentifier", "ortholog_gene.primaryIdentifier",
            "ortholog_organism.shortName", "relationship"
        )
        query.add_constraint("gene.primaryIdentifier", "CONTAINS", str(gName), code = "A")
        query.add_constraint("organism.shortName", "=", "A. thaliana", code = "B")
        query.add_constraint("ortholog_organism.shortName", "=", "E. salsugineum", code = "C")
        
        for row in query.rows():
            print(row["gene.primaryIdentifier"], row["ortholog_gene.primaryIdentifier"], 
                    row["ortholog_organism.shortName"], row["relationship"])
