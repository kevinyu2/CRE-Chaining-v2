import math

'''
Takes in a list of files <input_files> to get the ACRs. Randomly splits into two 
sets. Set 1 will contain <ref_fraction> of the ACRs.
'''
def generate_sets(input_files, output1_file, output2_file, ref_fraction):
    ACRs = []
    for file_path in input_files:
        with open(file_path) as file:
            for line in file:
                if line != "\n":
                    ACRs.append(line)
    ref_num = math.floor(len(ACRs) * ref_fraction)
    
    with open(output1_file, "w") as output:
        for i in range(ref_num):
            output.write(ACRs[i])
    
    with open(output2_file, "w") as output:
        for i in range(ref_num, len(ACRs)):
            output.write(ACRs[i])