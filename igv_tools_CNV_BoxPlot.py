#Coder: Vyomesh Javle
#Date: 07-June-22 to 11-June-22, 25-Oct-22 to 06-Nov-22

import os
from pathlib import Path
import time
import subprocess
import glob
import pandas as pd
import re
import plotly.express as px

print("\n Let's Get the Mutation Burden \n")

print(" ########## Taking Sample ID's from list_s.txt for Somatic DNA Samples  ########## ")
file1=open('list_s.txt', 'r')
data=file1.read()
sample_id_list_1=data.split("\n")
file1.close()

print(" ########## Taking Sample ID's from list_g.txt for Germline Samples  ########## ")
file2=open('list_g.txt', 'r')
data=file2.read()
sample_id_list_2=data.split("\n")
file2.close()

print(" ########## Taking Sample ID's from list_r.txt for Somatic RNA Samples  ########## ")
file3=open('list_r.txt', 'r')
data=file3.read()
sample_id_list_3=data.split("\n")
file3.close()

print(" ########## Taking Basespace Project Names and system name from list_project_names.txt  ########## ")
file4=open('list_project_names.txt', 'r')
data=file4.read()
sample_id_list_4=data.split("\n")
file4.close()
project_1=sample_id_list_4[0]
project_2=sample_id_list_4[1]
project_3=sample_id_list_4[2]
system_name=sample_id_list_4[3]
col_number=sample_id_list_4[4] 
ot_folder=sample_id_list_4[5]
print(" Project Somatic DNA : "+project_1+"\n"+" Project Germline    : "+project_2+"\n"+" Project Somatic RNA : "+project_3+"\n"+" Output folder name : "+ot_folder+"\n"+" System name : "+system_name+"\n")

local_path = os.getcwd()
mkdir_cmd_2='mkdir '+ local_path+"/"+ot_folder
os.system(mkdir_cmd_2)
loc_to_save=local_path+"/"+ot_folder

#enter the coloumns to be considered in final output file.
#Sample_ID,Chrom,Pos,Total_Count,Base_A,Base_G,Base_T,Base_C,Base_N,Base_INS,Base_DEL,Percent_Allelic_Burden_A,Percent_Allelic_Burden_G,Percent_Allelic_Burden_C,Percent_Allelic_Burden_T,Percent_Allelic_Burden_N,Percent_Allelic_Burden_INS,Percent_Allelic_Burden_DEL for input in nxt line 
file5=open('list_of_output_coloumns.txt', 'r')
data=file5.read()
Col_list=data.split("\n")
file5.close()

print("Enter the category number for column to be included in final file: 0 (Default)" )
col_num = int(col_number) #changed
col_ids=Col_list[col_num]
col_ids_list= col_ids.split(",")

print(" ########## Taking Mutation Location  ########## ")
Cap_Kit_data = pd.read_csv('/path/Region_file.csv', sep=',')

for i in range(len(Cap_Kit_data)):
    Kit_Reg_Num=Cap_Kit_data.loc[i, "Cap_Region_Num"]
    Kit_Reg_Cat=Cap_Kit_data.loc[i, "Cap_Region_Cat"]
    Gene = Cap_Kit_data.loc[i, "Gene_Name"] #need to check
    Start=Cap_Kit_data.loc[i, "Start"]
    End=Cap_Kit_data.loc[i, "End"]
    mut_loc=Cap_Kit_data.loc[i,"Combined"]
    range_mut_loc=mut_loc
    print("\n ######################## Starting the count & allelic Burden calculation for Somatic DNA samples "+range_mut_loc+"  ######################## \n")
    for ids_1 in sample_id_list_1:
        if ids_1 != "NA":
            print("Identifying the required mutation for: "+ids_1)
            comd_1="/home/"+system_name+"/Programs/IGV_2.13.0/igvtools count -w 1 --bases --query "+range_mut_loc+" /home/"+system_name+"/path/"+project_1+"/path/"+ids_1+"/Files/"+ids_1+".bam"+" "+loc_to_save+"/"+ids_1+"_"+Kit_Reg_Cat+"_"+str(Kit_Reg_Num)+"_"+range_mut_loc+"_.wig hg19"
            os.system(comd_1)
            print("Wig file created for: "+ids_1) 
        else:
            print("\n ######################## No ID Entered/Identified for Somatic DNA Samples ######################## \n")


    print("\n ######################## Starting the count & allelic Burden calculation for Germline samples "+range_mut_loc+" ######################## \n")
    for ids_2 in sample_id_list_2:
        if ids_2 != "NA":
            print("Identifying the required mutation for: "+ids_2) 
            comd_2="/home/"+system_name+"/Programs/IGV_2.13.0/igvtools count -w 1 --bases --query "+range_mut_loc+" /home/"+system_name+"/path/"+project_2+"/path/"+ids_2+"/Files/"+ids_2+".bam"+" "+loc_to_save+"/"+ids_2+"_"+Kit_Reg_Cat+"_"+str(Kit_Reg_Num)+"_"+range_mut_loc+"_.wig hg19" 
            os.system(comd_2) 
            print("Wig file created for: "+ids_2)
        else:
            print("\n ######################## No ID Entered/Identified for Germline Samples ########################  \n")

    print("\n ######################## Starting the count & allelic Burden calculation for Somatic RNA samples "+range_mut_loc+"  ########################  \n")
    for ids_3 in sample_id_list_3:
        if ids_3 != "NA": 
            print("Identifying the required mutation for: "+ids_3) 
            comd_3="/home/"+system_name+"/Programs/IGV_2.13.0/igvtools count -w 1 --bases --query "+range_mut_loc+" /home/"+system_name+"/path/"+project_3+"/path/"+ids_3+"_*/Files/"+ids_3+".bam"+" "+loc_to_save+"/"+ids_3+"_"+Kit_Reg_Cat+"_"+str(Kit_Reg_Num)+"_"+range_mut_loc+"_.wig hg19"
            os.system(comd_3) 
            print("Wig file created for: "+ids_3) 
        else:
            print(" \n ######################## No ID Entered/Identified for Somatic RNA Samples ########################  \n")
    print(" \n #################### All Samples with "+range_mut_loc+" Done ####################  \n")

print(" \n #################### Initiating the compilation of CNV region Total Count in the Final Excel file for "+range_mut_loc+"  ####################  \n")

###### All WIG files taken together for CNV's ########
files = glob.glob(os.path.join(loc_to_save,"*.wig"))
print(files)

files_1=[]
files_zero=[]
for fp2 in files:
    if os.path.getsize(fp2) > 0:
        files_1.append(fp2)
    else:
        print("\n Temporary Skipping the file(s) "+os.path.basename(fp2)+" as files are empty.\n")
        files_zero.append(fp2)
     
df_n2=pd.DataFrame(columns=['Pos','Base_A','Base_C','Base_G','Base_T','Base_N','Base_INS','Base_DEL','Chrom','Sample_ID'],dtype=object)

proceed_1="YES"
print("Proceed to modify and append the wig Files & Calculating Allelic Burden: YES (Default)" ) 
if proceed_1.lower() == "yes":
    for fp in files_1:
        df= pd.read_table(fp, skiprows = [0,1,2],header=None,index_col=False)
        df_n1=df.rename(columns={df.columns[0]: 'Pos', df.columns[1]: 'Base_A', df.columns[2]: 'Base_C', df.columns[3]: 'Base_G', df.columns[4]: 'Base_T',df.columns[5]: 'Base_N',df.columns[6]: 'Base_DEL',df.columns[7]: 'Base_INS'})
        df_n1['Pos']=df_n1['Pos'].astype(str)
        range_mut_loc=str(os.path.basename(fp).split('_')[3])
        df_n1['Chrom']=range_mut_loc.split(':')[0] #work here if there is only one chromshome (need to change if Exon_data has multiple chromosomes) #VJ #26Aug22
        df_n1['Sample_ID']=os.path.basename(fp).split('_')[0]
        df_n1['Kit_Reg_Cat']=os.path.basename(fp).split('_')[1]
        df_n1['Kit_Reg_Num']=os.path.basename(fp).split('_')[2]
        df_n1['Depth']= df_n1['Base_A']+df_n1['Base_C']+df_n1['Base_G']+df_n1['Base_T']+df_n1['Base_N']+df_n1['Base_INS']+df_n1['Base_DEL']  
        range_mut_loc=str(os.path.basename(fp).split('_')[3])
        #if re.search(r'chr\d*:\d*-',range_mut_loc):
        range_1=range_mut_loc.split(':')[1]
        df_n2=df_n2.append(df_n1.loc[df_n1['Pos'].between(range_1.split('-')[0],range_1.split('-')[1])],ignore_index = True)
    df_n3=df_n2[col_ids_list]

print("CNV Total Count sheets is ready")

####### Creating the Box Plot w.r.t added 26-Sep-22 VJ ####
depth_file=loc_to_save+"/"+ot_folder+"_output.csv"  
ot_total_count = pd.read_csv(depth_file, sep=',',usecols=['Sample_ID', 'Chrom', 'Pos', 'Depth', 'Kit_Reg_Num','Kit_Reg_Cat'])
rslt_ot_total_count  = ot_total_count.sort_values(by = ['Kit_Reg_Num', 'Sample_ID'])
rslt_ot_total_count_2 = rslt_ot_total_count[rslt_ot_total_count.Depth >= 100]

#For Multiple samples
fig = px.box(rslt_ot_total_count_2, x ="Kit_Reg_Cat", y="Depth", notched=True, color="Sample_ID",title="CNA Amplification Plot for SULF2") #need to change
a=loc_to_save+"/"+ot_folder+"_Amplification1.html"
fig.write_html(a)

######

