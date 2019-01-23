# Dictionary with strings to replace and what to replace them with
import sys
import os
import tempfile
import glob
import re
import io
from Bio import SeqIO



replace_strings = {}
with open("Sequencing_strain_PM.txt", "r") as id_file:
    # Read file line-by-line
    for line in id_file.readlines():
        # Split line on TAB
        ids = line.strip().split("\t")
        # Fist entry is the original ID
        original_id = ids[0]
        # Second entry is your ID
        my_id = ids[1]
        # Add both to our dictionary of strings to replace
        replace_strings[original_id] = my_id
        #print(replace_strings)

        #alternation code
    #for line in id_file:
     #   ids = line.strip().split('\t')
      #  replace_strings[ids[0]] = ids[1]

# Read file with sequences, loop over all gbk files
for file in glob.iglob('*.gbk'):
    with open(file, "r+") as infile:
    # Read each line of file into a list
        content = infile.readlines()
    # Keep a list of the lines with the replaced strings
        new_content = []
    # Loop lines in the file content
        for line in content:
            new_line = line
        # Find and replace any original_id with your own ids in the line of content and add it to our list of replaced lines
            for original_id, my_id in replace_strings.items():
                new_line = new_line.replace(original_id, my_id)
            new_content.append(new_line)
    a=file.rsplit('.', 1)[0]
    # Write replaced content to a new file called "outfile.txt"
    with open(a+"_pm"+".gbk", "w") as outfile:
        for line in new_content:
            outfile.write(line)



# i in xrange(10):
 #   with io.open("file_" + str(i) + ".dat", 'w', encoding='utf-8') as f:
  #     f.write(str(func(i))