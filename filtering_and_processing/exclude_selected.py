#import re
import string
import sys
#import gzip

def write_filtered_file(file_name,to_write):
    with open(file_name, "w") as file_handle:
        file_handle.write("label,pixelarea,probecount,expid\n")
        for line in to_write:
            file_handle.write(line+'\n')

def read_filtering_line(file_name):
    "tab delimited file. must be sorted by image"
    for line in open(file_name, "r"):
        fields=line.split('\t')
        if fields[0].strip()!="Image":
            image=int(fields[0].strip())
            exclude_str=fields[3].strip('\"')
            if exclude_str=="":
                print("empty exclusion for image "+str(image))
                exclude=[]
            else:
                exclude_fields=exclude_str.split(',')
                exclude=[int(x.strip()) for x in exclude_fields]
            stuck_px_str=fields[4].strip('\"')
            if stuck_px_str=="":
                print("empty stuck_px for image "+str(image))
                stuck_px=[]
            else:
                stuck_px_fields=stuck_px_str.split(',')
                stuck_px=[int(x.strip()) for x in stuck_px_fields]
            yield (image,exclude,stuck_px)

def read_result_line(file_name):
    "csv file output of FISHTOOLS. must be sorted by expid and then by label"
    for line in open(file_name, "r"):
        fields=line.split(',')
        if fields[0].strip()!="label":
            label=int(fields[0].strip())
            pixelarea=int(fields[1].strip())
            probecount=int(fields[2].strip())
            expid=int(fields[3].strip())
            yield (label,pixelarea,probecount,expid)


def main(file_in,file_out):
    to_write=[]
    filename_out=file_out.split(".csv")
    filtered_file_name=filename_out[0]+"_filtered.csv"
    filter_line=read_filtering_line(file_in)
    image,exclude,stuck_px=next(filter_line)
    for label,pixelarea,probecount,expid in read_result_line(file_out):
        if expid!=image:
            image,exclude,stuck_px=next(filter_line)
        if expid!=image:
            print("Error! Image "+str(image)+" is not "+str(expid))
        if label in stuck_px:
            probecount-=1
        if label in exclude:
            print("removed "+str(label)+" from image "+str(expid))
        else:
            to_write.append(str(label)+','+str(pixelarea)+','+str(probecount)+','+str(expid))
    write_filtered_file(filtered_file_name,to_write)



main(sys.argv[1],sys.argv[2])