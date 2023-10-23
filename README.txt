README for FragAssemble.java

Author: Nick Fisk

Usage: FragAssemble {name of FASTA Formated sequence file}{minimum_overlap}

IMPORTANT!!!!!!!
Each FASTA entry in the FASTA file must be uniquely named
However, duplicate entries are handled fine.
If two entries have the same name and are not duplicate sequences,
the program will NOT assemble the fragments in the described way.


FragAssemble relies on Contig.java and Frag.java to operate.
It assembles DNA fragments based on a greedy beam search. 
It condenses fragments into single fragments (while keeping
the appropriate inter-fragment relations) as the graph is
being traversed, so the runtime/operation should decrease as
the program runs on. The contigs are numbered arbitrarily
by the order they are added to the list of contigs. 

The results are printed to the standard out. There is no
prompting for user input, so it is safe to redirect the output
to a file, if desired.  
