from fasta_to_physicochemicalprops_ import FastaToPhysicochemical3mers


input_file_path = "input_file.fasta"
folder_for_output = "./output/"
encoder = FastaToPhysicochemical3mers(folder_for_output=folder_for_output)
encoder.convert_to_physicochemical_properties(input_file_path)
encoder.convert_to_physicochemical_stats(input_file_path)