1.1.3
=====

* Small BUG FIX in the progress bar that caused epochs to be misrepresented

1.1.2
=====

* Genome import is now much more memory efficient

* BUG FIX: find in BinarySequence could not find subsequences at the tail of sequence

* Added a built-in datawrap with only chr1 and y

* Readme update with more infos about importation and link to doc
 
1.1.1
=====

Much better SNP filtering interface
------------------------------------
Easier and much morr coherent:

* SNP filtering has now its own module

* SNP Filters are now objects

* SNP Filters must return SequenceSNP, SNPInsert, SNPDeletion or None objects

1.0.0
=====
Freshly hatched

