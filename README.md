# Package need to be downloaded
(under environment Python3.8)
- biopython 1.78
- gffutils 0.10.1
- pandas 1.1.4
- numpy 1.19.4
- pip 20.2.4

# Q1a
put `Vibrio_cholerae.GFC_11.dna.toplevel.fa` , `Vibrio_cholerae.GFC_11.37.gff3` and `Q1a.py` in the same working directory

# Q1b
Make sure `config.py`, `Q1b.py` and`Vibrio_cholerae.GFC_11.dna.toplevel.fa` are on the same working directory<br/>
To run the two input files, on the configuration settings of `Q1b.py`, find Parameters and entre `"Vibrio_cholerae.GFC_11.dna.toplevel.fa" "config.py"`. 

# Q1c
Make sure `config.py`, `Q1b.py` and`Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa` are on the same working directory<br/>
On the configuration settings of `Q1b.py`, find Parameters and entre `"Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa" "config.py"`<br/>
To get output name changed, change code on `Q1b.py` line 27 to `output = open(r"1c_out.gff3", "w")`

# Q1d
Make sure `Q1d.py`,`Vibrio_vulnificus.ASM74310v1.37.gff3` and `1c_out.gff3` are on the same working directory. 


