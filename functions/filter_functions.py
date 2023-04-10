""" 
Script to filter mmseqs easy_taxonomy functions 

Selects the best unirefID for each read below some e value cutoff and ignoring the uniparc ids 
""" 

# imports 
import click 
import pickle 
import csv 

# click 
@click.command()
@click.option('--functions_path', '-fp', type=click.Path(exists=True), required = True, help='path to functions output from mmseqs')
@click.option('--out_path', '-o', type=click.Path(exists=False), required = True, help='path to the output dictionary')
@click.option('--e_cutoff', '-e', type=float, help = 'e value cutoff')

def main(functions_path, out_path, e_cutoff):
    """
    Filter the mmseqs filter output 
    """
    
    # dictionary of functions 
    data = {} 
    
    # variables to keep contol of the data on hand 
    current_read = '' 
    best_evalue = '' 
    best_id = ''

    with open(functions_path, mode='r') as file:
    
        reader = csv.reader(file, delimiter='\t')

        for row in reader:

            #if the read has an e-value above a threshold and is not from UniParc then use it 
            # in this section loop through each row 

            read = row[0]
            e_val = float(row[-2]) 
            this_id = row[1]

            #only consider if it is not a uniparc id 
            if this_id[:11] != 'UniRef50_UP':

                #if a new read is encountered 
                if current_read != read:

                    #update for the new read 
                    current_read = read 
                    best_evalue = e_val

                    if e_val <= e_cutoff: 
                        data[read] = this_id

                #if the read has been seen previously 
                else:            

                    #update if the e_val is better than the current best and not a uniparc id 
                    if e_val < best_evalue:

                        best_evalue = e_val

                        if e_val <= e_cutoff: 
                            data[read] = this_id

    # save the dictionary to a file 
    with open(out_path, 'wb') as handle: 
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
if __name__ == '__main__':
    main()
