git clone https://github.com/tseemann/barrnap.git
cd barrnap/bin
./barrnap --help

# Rename all *.txt to *.text
for f in *.fasta; do 
    mv -- "$f" "${f%.fasta}.fna"
done
#loop all files in folder$.fasta
#!/bin/bash
for filename in `ls *.fna`; 
   do
     for i in {0..4}; 

       do fold -w70 $filename>a$filename;
        
     done 
 done

 awk '/>/{print ">Cdifficile Prochlorococcus CIM-523" ++i; next}{print}' < aCIM-522.fna>bCIM-522.fna


for filename in `ls *.fna`; 

       do barrnap --outseq out$filename $filename ; 
 done

barrnap --outseq 11.fasta CIM-072.fna