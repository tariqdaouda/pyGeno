Querying
=========

pyGeno is a personal database that you can query un many ways. Special emphasis has been placed upon ease of use, and you only need to remember two functions::

	* get()
	* help()

**get()** can be called from any pyGeno object to get any other object 

**help()** is you best friend when you get lost using **get()**. When called, it will give the list of all field that you can use in get queries. You can call it either of the class::

	Gene.help()

Or on the object::

	ref = Genome(name = "GRCh37.75")
	g = ref.get(Gene, name = "TPST2")[0]

	g.help()

Both will print::

	'Available fields for Gene: end, name, chromosome, start, biotype, id, strand, genome'
