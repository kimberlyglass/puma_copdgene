
import numpy as np
from linecache import getline, clearcache
import os
import pandas as pd
import os



path_files_by_subj = snakemake.params.path_puma_sbj
path_files_by_mirna = snakemake.params.path_puma_mirna # path_puma_files + 'by_mirna/'
mirna = snakemake.wildcards.mirna# #snakemake.wildcards.mirna


tmp_list = []
print(mirna)
def get_mirna_k(path_file,mirna):
    # Read only the first column
    tmp_list_mirna = pd.read_csv(path_file, usecols=[0])
    k_th_mirna = np.where(np.array(tmp_list_mirna)==mirna)[0][0] 
    return k_th_mirna

list_files = os.listdir(path_files_by_subj)[::-1]


for  i,file_name in enumerate(list_files):
    sbj = file_name.split('.')[0]
    if i==0:
        k_th_mirna = get_mirna_k(path_files_by_subj+file_name,mirna)
        columns_name0 = getline(path_files_by_subj+file_name,1)
        tmp_list = [columns_name0]
        line = getline(path_files_by_subj+file_name,k_th_mirna+2)
        try:
            os.remove(path_files_by_mirna+'%s.csv'%mirna)
        except:
            continue
    else:
        columns_name = getline(path_files_by_subj+file_name,1)
        assert columns_name0 == columns_name, 'wrong columns order'
        line = getline(path_files_by_subj+file_name,k_th_mirna+2)
    mirna_line = line.split(',')[0]
    if mirna_line!=mirna:
        print('Updated k-th')
        # mirna is different recompute
        k_th_mirna = get_mirna_k(path_files_by_subj+file_name,mirna)
        line = getline(path_files_by_subj+file_name,k_th_mirna+2)
        mirna_line = line.split(',')[0]
    # make sure new match
    assert mirna_line==mirna,'mirna difference '+mirna_line+' '+mirna 

    line = line.replace(mirna,sbj,1) # replace only first occurrence
    tmp_list.append(line)
    clearcache()
    if i%10==0:#rate is not None and rate < rate0/2:
        print(i,'Save file')

        #if list getting to big, save it
        with open(path_files_by_mirna+'%s.csv'%mirna,'a') as f:
            f.write(''.join(tmp_list))
        
        tmp_list = []
        
        
    """files=( {path_input_lio_file}) #/d/tmp/remge/PUMA_final/Lioness_Puma/by_subj/*.csv )
        total_files=${#files[@]}

        for i in "${!files[@]}"; do
            # Process the file using awk
                                    # Remove the .csv extension
            if [[ $i == 0 ]]; then
            # Get the first element of the first file
            header = $(awk -F ',' 'NR==1 {print $0}' "${files[$i]}")
            list_mirna=($(awk -F ',' 'NR>1 {print $1}' "${files[$i]}"))
            # Iterate over the array and create a file for each element
            for mirna in "${list_mirna[@]}"; do
                echo $header > "/d/tmp/remge/PUMA_final/Lioness_Puma/by_mirna/${mirna}.txt"
            done

            # Copy the first line of the first file to all files named with the first element
            # This assumes the first line of the first file contains data like "a 3 1 5"
            
            fi  
            awk -F',' '{
            # Extract the file name without the extension
            split(FILENAME, path_parts, "/");           # Split the path by '/'
            file = path_parts[length(path_parts)];      # Get the last part (file name with extension)
            sub(".csv$", "", file);                     # Remove the .csv extension

            # Print the extracted content to the corresponding output file
            print file "," substr($0, index($0,$2)) >> "/d/tmp/remge/PUMA_final/Lioness_Puma/by_mirna/" $1 ".csv";
            }' "${files[$i]}"
            
            # Print progress
            echo "Processed file $((i+1)) of $total_files: ${files[$i]}"
        done"""