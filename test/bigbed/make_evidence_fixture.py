#!/usr/bin/env python3
"""Regenerate samples/human.chr14.longintron.introns.bb (needs pybigtools).
Two Intron evidence records matching the human_intron fixture slice, used by the
`human_intron_bb` regression case to exercise the bigBed -R evidence path."""
import pybigtools
b = pybigtools.open("samples/human.chr14.longintron.introns.bb", "w")
b.write({"human.chr14.longintron": 250001},
        iter([("human.chr14.longintron", 41750, 70730, ".\t100\t-\tIntron"),
              ("human.chr14.longintron", 70789, 119365, ".\t100\t-\tIntron")]))
b.close()
print("wrote samples/human.chr14.longintron.introns.bb")
