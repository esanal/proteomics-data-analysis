####    Functions for data editing    ####


##Functions in module#

#annotate(path,column,output)
#filter_by_column(path,column_name,what_to_filter_out,output_name)
#normalize_by_column(path,column_name,normal_info_column,log_type,output_name)



####---------------------------------####


####Codes####


####Go annotation####
def annotate(path,column,output):
    from bioservices import QuickGO

    s = QuickGO()
    
    data = open(path,'r')
    data_o = open(output,'w')
    
    header = data.readline()
    new_header = header.strip('\n')+'\tMembrane\tCytoplasm\tOthers\tAnnotation\n'
    data_o.write(new_header)
    
    #find the index of column to be filtered   
    i = header.strip('\n').split('\t').index(column)
    
    for j,line in enumerate(data):
        print j," element is done         \r",
        line_s = line.split('\t')
        p_id = line_s[i]
        
        ann = s.Annotation_from_protein(protein = p_id, frmt = 'tsv', tax=9606, source = 'UniProt', col = 'goName')
        ann_list = list(set(ann.goName))
        pm = 'FALSE'
        cyt = 'FALSE'
        others = 'TRUE'
    
        if any('membrane' in s for s in ann_list):
            pm = 'TRUE'
            others = 'FALSE'
        elif 'cytosol' in ann_list or 'cytoplasm' in ann_list:
            cyt = 'TRUE'
            others='FALSE'
    
        
        line_w = line.strip('\n')+'\t'+pm+'\t'+cyt+'\t'+others+'\t'+','.join(ann_list)+'\n'    
        data_o.write(line_w)
        
    data.close()
    data_o.close()
    
    
    

####Filtering by specified column####
def filter_by_column(path,column_name,what_to_filter_out,output_name):
    
    #open data file to filter and file to write
    data = open(path,'r')
    output = open(output_name,'w')
    
    #get the header
    header = data.readline()
    output.write(header)
    #find the index of column to be filtered   
    i = header.strip('\n').split('\t').index(column_name)
    
    
    #check everyline if what_to_filter is in there
    for j,line in enumerate(data):
        print j," element is done         \r",
        if what_to_filter_out == line.strip('\n').split('\t')[i]:
            continue
        else:
            output.write(line)
        
    #close files 
    data.close()
    output.close()
    
    
    

####Normalize data####
def normalize_by_column(path,column_name,normal_info_column,log_type,output_name):
    
    import numpy as np
    import math
    #open data file to norm. and file to write
    data = open(path,'r')
    output = open(output_name,'w')
    
    #get the header
    header = data.readline()
    output.write(header.strip('\n')+'\t'+str(log_type)+' Normalized '+column_name+'\n')
    
    #find the index of column to be normalized   
    i_values = header.strip('\n').split('\t').index(column_name)
        
    #find the index of column with normal info  
    i_inf = header.strip('\n').split('\t').index(normal_info_column)

    #array for normal data
    normal_data = []
    
    #loop to get normal data 
    for line in data:
        line_s = line.strip('\n').split('\t')
        if line_s[i_inf] == "TRUE":
            continue
        else:
            normal_data.append(float(line_s[i_values]))
            
    #find denominator for normalization
    denominator = np.median(normal_data)
    
    #loop to normalize
    for j,line in enumerate(data):
        print j," element is done         \r",
        line_s = line.strip('\n').split('\t')
        norm_val = math.log(float(line_s[i_values]) / denominator,log_type)
        output.write(line.strip('\n')+'\t'+str(norm_val)+'\n')
        
    #close files 
    data.close()
    output.close()




####Create "Uniprot ID" column####
def make_uniprot(path, col_name, output):
    
    data = open(path,'r')
    output = open(output,'w')
    
    #get the header
    header = data.readline()
    output.write(header.strip('\n')+'\t'+'Uniprot ID'+'\n')
    
    #find the index of column to be normalized   
    i = header.strip('\n').split('\t').index(col_name)
    
    #loop
    for j,line in enumerate(data):
        print j," element is done         \r",
        line_s = line.strip('\n').split('\t')
        string = line_s[i]
        f_l = string.find('|')
        s_l = string[f_l+1:].find('|') + f_l
        u_id = string[f_l+1:s_l+1]
        
        output.write(line.strip('\n')+'\t'+u_id+'\n')
        
        
    
        

    
            