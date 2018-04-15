help_description = ""
flag = 0 
for line in open("README.md"):
    if flag == 1: 
        if "````" in line: 
            pass 
        else:
            help_description += line 

    if "````help" in line:
        flag = 1
