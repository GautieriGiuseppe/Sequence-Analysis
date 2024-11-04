Multiple Sequence Alignment script
It aligns two or more sequences using two of the most used tools to perform MSA, which are ClustalW and Muscle.
Starting from the list of the paths in which FASTA sequences are stored, it combines all the sequences in a single file and then runs clustalw or muscle using conda and subprocess module.
In order to make this work you need to have conda installed because it makes usage of the terminal.
You can switch from one tool to the other simply changing the string passed to the function 'perform_msa', the default one is clustalw.
