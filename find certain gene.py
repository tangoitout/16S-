# Dictionary with strings to replace and what to replace them with
replace_strings = {}
with open("replace_ids.tsv", "r") as id_file:
    # Read file line-by-line


	for line in id_file:
    	ids = line.strip().split('\t')
    	replace_strings[ids[0]] = ids[1]



# Read file with sequences, called "sequence.txt"
with open("sequence.txt", "r+") as infile:
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

    # Write replaced content to a new file called "outfile.txt"
    with open("outfile.txt", "w") as outfile:
        for line in new_content:
            outfile.write(line)


