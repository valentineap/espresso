# Command to create a new example folder:
# In ROOT, execute:
# python utils/new_contribution/create_new_contrib.py <example-name>
# Replacing <example_name> with the new example name.

from glob import glob
import sys
import os
from shutil import copyfile
from pathlib import Path


from click import UsageError

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles


current_dir = Path(__file__).parent
root_dir = current_dir.parent.parent
CONTRIB_FOLDER = str(root_dir / "contrib")
TEMPLATE_FOLDER = str(root_dir / "utils/new_contribution/_template")

def main():
    print("🥰 Thanks for contributing! \nWe are generating new contribution component from template...\n")

    # validate example name
    if len(sys.argv)!=2:
        raise RuntimeError("No example name detected.\n\nUsage: python create_new_contrib.py EXAMPLE_NAME\n\n")
    example_name = sys.argv[-1]
    existing_examples = glob(CONTRIB_FOLDER+"/*/")
    existing_examples = [e for e in existing_examples]
    existing_examples = [e.split("/")[1] for e in existing_examples]
    if example_name in existing_examples:
        raise ValueError("The example name provided already exists, please choose another name")

    # convert example name to other formats
    example_name_capitalised = example_name.title().replace("_", " ")
    example_name_no_space = example_name_capitalised.replace(" ", "")

    # make new folders and subfolders
    new_subfolder = CONTRIB_FOLDER+"/"+example_name
    os.makedirs(new_subfolder)    
    os.makedirs(new_subfolder+"/data")    

    # generate file names
    template_files=getListOfFiles(TEMPLATE_FOLDER)
    new_files = [f.replace(TEMPLATE_FOLDER, new_subfolder) for f in template_files]
    new_files = [f.replace("example_name", example_name) for f in new_files]

    # generate example subfolder from template
    print("🗂  Copying files...")
    for template_file, new_file in zip(template_files, new_files):
        print("file: "+template_file+" -> "+new_file)
        if "README" in template_file:
            with open(template_file, "r") as template_f:
                content = template_f.read()
            content = content.replace("example_name", example_name)
            content = content.replace("Example Name Title", example_name_capitalised)
            with open(new_file, "w") as new_f:
                new_f.write(content)
        elif "__init__.py" in template_file:
            with open(template_file, "r") as template_f:
                content = template_f.read()
            content = content.replace("example_name", example_name)
            content = content.replace("ExampleName", example_name_no_space)
            with open(new_file, "w") as new_f:
                new_f.write(content)
        elif "example_name.py" in template_file:
            with open(template_file, "r") as template_f:
                content = template_f.read()
            content = content.replace("ExampleName", example_name_no_space)
            with open(new_file, "w") as new_f:
                new_f.write(content)
        else:
            copyfile(template_file, new_file)
    print("\n🎉 OK. Please navigate to " + new_subfolder + " to write your own example. ")

if __name__ == '__main__':
    main()
