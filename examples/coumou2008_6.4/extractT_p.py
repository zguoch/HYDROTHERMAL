import numpy as np 
import linecache
import matplotlib 

TAG_COORD_X='X-Direction Node Coordinates'
TAG_COORD_Y='Y-Direction Node Coordinates'
TAG_TIME='Simulation time:'
TAG_TEMPERATURE='--- Temperature Values ---'
TAG_COMPLETED='Simulation Completed'
def extract_field(ALLDATA,start):
    i=start
    field=[]
    while(not (TAG_TIME in ALLDATA[i])):
        i=i+1
        if(TAG_COMPLETED in ALLDATA[i]):
            break
        str_array=np.array(ALLDATA[i].replace('\n','').split(' '))
        ind_num=(str_array!='')
        # print(str_array)
        coord_array=np.array(str_array[ind_num],dtype=float)
        field.extend(coord_array[1:])
        i=i+4
    # print(i,field)
    return field,i
def extract_T_1D(fname):
    result_ht={}
    result_ht['t']=[]
    result_ht['field']=[]
    result_ht['x']=[]
    result_ht['z']=[]

    ALLDATA=linecache.getlines(fname)
    for i in range(0,len(ALLDATA)):
        str_line=ALLDATA[i]
        # 1. x coordinate
        if(TAG_COORD_X in str_line):
            result_ht['unit_x']=str_line.split('(')[-1].split(')')[0]
            i=i+2
            while(ALLDATA[i]!='\n'):
                str_array=np.array(ALLDATA[i+1].replace('\n','').split(' '))
                ind_num=(str_array!='')
                coord_array=np.array(str_array[ind_num],dtype=float)
                result_ht['x'].extend(coord_array)
                i=i+5
        # 2. z coordinate
        if(TAG_COORD_Y in str_line):
            result_ht['unit_y']=str_line.split('(')[-1].split(')')[0]
            i=i+2
            while(ALLDATA[i]!='\n'):
                str_array=np.array(ALLDATA[i+1].replace('\n','').split(' '))
                ind_num=(str_array!='')
                coord_array=np.array(str_array[ind_num],dtype=float)
                result_ht['z'].extend(coord_array)
                i=i+5
        # 3. extract field for every write time
        if(TAG_TIME in str_line):
            t=np.float(str_line.split(TAG_TIME)[-1].split('(')[0])
            unit_t=str_line.split(TAG_TIME)[-1].split('(')[-1].split(')')[0]
            result_ht['unit_t']=unit_t
            result_ht['t'].append(t)
            i=i+2
            str_line=ALLDATA[i]
            result_ht['unit_T']=str_line.split('(')[-1].split(')')[0]
            field,end=extract_field(ALLDATA,start=i+2)
            result_ht['field'].append(field)
            i=end
            # print(field)
    linecache.clearcache()
    if(result_ht['t']!=[]):
        result_ht['t']=np.array(result_ht['t'],dtype=float)
    if(result_ht['field']!=[]):
        result_ht['field']=np.array(result_ht['field'],dtype=float)
    if(result_ht['x']!=[]):
        result_ht['x']=np.array(result_ht['x'],dtype=float)
    if(result_ht['z']!=[]):
        result_ht['z']=np.array(result_ht['z'],dtype=float)
    return result_ht

def write2file(TData,pData):
    x=TData['x']
    for t,p,T in zip(TData['t'],pData['field'],TData['field']):
        fpout=open('T_p_'+str(t)+'.txt','w')
        fpout.write('# x\tT\tp\n')
        for x0, T0, p0 in zip(x,T,p):
            fpout.write('%.6f\t%.6f\t%.6f\n'% (x0,T0,p0))
        fpout.close()


# ------------extract data-------------------
fname_T='Out_temperature.dx20'
fname_p='Out_pressure.dx20'
TData=extract_T_1D(fname_T)
pData=extract_T_1D(fname_p)
write2file(TData,pData)

