import math
import scipy
import time
import numpy as np

def save_shakemap_to_xml(eventID,IMTs,exposureLocations,data):

    validIMTs = ['PGA', 'SA(0.3)', 'SA(1.0)', 'SA(3.0)']
    noLocations = len(exposureLocations)
    gmfs = np.zeros((noLocations*len(IMTs),1))
    data = np.array(data)
    for i in range(4):
        if validIMTs[i] in IMTs:
            for j in range(len(data)):
                gmfs[i*noLocations+j][0]=math.exp(data[j,2+i])

    save_gmfs_to_xml(eventID,IMTs,exposureLocations,gmfs)

def save_gmfs_to_xml(eventID,IMTs,exposureLocations,gmfs):

    noLocations = len(exposureLocations)
    validIMTs = ['PGA', 'SA(0.3)', 'SA(1.0)', 'SA(3.0)']

    xml = open(str(eventID)+'.xml',"w")
    xml.write('<?xml version="1.0" encoding="utf-8"?>\n')
    xml.write('<nrml\n')
    xml.write('xmlns="http://openquake.org/xmlns/nrml/0.4"\n')
    xml.write('xmlns:gml="http://www.opengis.net/gml"\n')
    xml.write('>\n')
    xml.write(' <gmfCollection\n')
    xml.write(' gsimTreePath=""\n')
    xml.write(' sourceModelTreePath=""\n')
    xml.write(' >\n')
    xml.write('         <gmfSet\n')
    xml.write('         stochasticEventSetId="1"\n')
    xml.write('         >\n')

    for i in range(4):
        if validIMTs[i] in IMTs:
            for j in range(len(gmfs[0])):
                xml.write('                     <gmf\n')
                if validIMTs[i] == 'PGA':
                    xml.write('                 IMT="PGA"\n')
                elif validIMTs[i][0:2] == 'SA':
                    T = float(validIMTs[i].replace("SA(","").replace(")",""))
                    xml.write('                 IMT="SA"\n')
                    xml.write('                 saDamping="5.0"\n')
                    xml.write('                 saPeriod="'+str(T)+'"\n')
                xml.write('                     ruptureId="scenario-'+str(j)+'"\n')
                xml.write('                     >\n')

                for iloc in range(noLocations):
                    xml.write('                         <node gmv="'+str(gmfs[i*noLocations+iloc][j])+'" lat="'+str(exposureLocations[iloc][1])+'" lon="'+str(exposureLocations[iloc][0])+'"/>\n')
                xml.write('                     </gmf>\n')
    xml.write('         </gmfSet>\n')
    xml.write(' </gmfCollection>\n')
    xml.write('</nrml>')
    xml.close()

def calculate_percentiles_gmfs(percentiles,IMTs,exposureLocations,gmfs):

    noLocations = int(len(gmfs)/len(IMTs))
    perGMFs = np.zeros((noLocations*len(IMTs),1))

    for iper in range(len(percentiles)):
        for i in range(len(IMTs)):
            for iloc in range(noLocations):
                perGMFs[i*noLocations+iloc]=np.percentile(gmfs[i*noLocations+iloc],percentiles[iper])
        save_gmfs_to_xml('gmf_'+str(percentiles[iper]),IMTs,exposureLocations,perGMFs)

def save_gmfs_to_csv(eventID,exposureLocations,gmfs):

    out_PGA = open(str(eventID)+'_PGA.csv',"w")
    out_Sa03 = open(str(eventID)+'_Sa03.csv',"w")
    out_Sa10 = open(str(eventID)+'_Sa10.csv',"w")
    out_Sa30 = open(str(eventID)+'_Sa30.csv',"w")
    noLocations = int(len(gmfs)/4)
    noGmvs = len(gmfs[0])

    out_PGA.write(buildHeader(noGmvs)+'\n')
    out_Sa03.write(buildHeader(noGmvs)+'\n')
    out_Sa10.write(buildHeader(noGmvs)+'\n')
    out_Sa30.write(buildHeader(noGmvs)+'\n')

    for iloc in range(noLocations):
        out_PGA.write(str(exposureLocations[iloc][0])+','+str(exposureLocations[iloc][1])+buildRow(gmfs[iloc])+'\n')
        out_Sa03.write(str(exposureLocations[iloc][0])+','+str(exposureLocations[iloc][1])+buildRow(gmfs[iloc+noLocations])+'\n')
        out_Sa10.write(str(exposureLocations[iloc][0])+','+str(exposureLocations[iloc][1])+buildRow(gmfs[iloc+2*noLocations])+'\n')
        out_Sa30.write(str(exposureLocations[iloc][0])+','+str(exposureLocations[iloc][1])+buildRow(gmfs[iloc+3*noLocations])+'\n')

    out_PGA.close()
    out_Sa03.close()
    out_Sa10.close()
    out_Sa30.close()

def buildHeader(noValues):

    row = 'lon,lat'
    for value in range(noValues):
        row = row +',gmf'+str(value)

    return row


def buildRow(gmvs):

    row = ''
    for gmv in gmvs:
        row = row +','+str(gmv)

    return row

def parse_nrml_gmf(gmfFile):

    IMTs = []
    locations = []
    gmvs = []

    file=open(gmfFile)
    lines=file.readlines()
    file.close

    for i in range(len(lines)):
        if lines[i].find('PGA')>0 or lines[i].find('SA')>0 :
            line = lines[i].strip('\n').strip(' ').strip('\t').split('"')
            if line[1] == 'PGA':
                IMT = line[1]
                firstLine = i+3
            if line[1] == 'SA':
                line2 = lines[i+2].strip('\n').strip('\t').split('"')
                IMT = line[1]+'('+line2[1]+')'
                firstLine = i+5
            line = lines[firstLine].strip('\n').strip(' ').strip('\t').split('"')
            counter = 0
            while line[0] == '<node gmv=':
                IMTs.append(IMT)
                locations.append([float(line[5]),float(line[3])])
                gmvs.append(float(line[1]))
                counter = counter+1
                line = lines[firstLine+counter].strip('\n').strip(' ').strip('\t').split('"')

    uniqueIMT = []
    uniquelocations = []
    for IMT in IMTs:
        if IMT not in uniqueIMT:
            uniqueIMT.append(IMT)
    for location in locations:
        if location not in uniquelocations:
            uniquelocations.append(location)

    noIMT = len(uniqueIMT)
    noLocations = len(uniquelocations)

    setGMFs = np.zeros((noIMT*noLocations,len(IMTs)/(noIMT*noLocations)))
    counterGMFs = np.zeros((noIMT*noLocations,1))

    for igmv in range(len(gmvs)):
        IMTId = uniqueIMT.index(IMTs[igmv])
        locId = uniquelocations.index(locations[igmv])
        gmfID = IMTId*noLocations+locId
        setGMFs[gmfID][counterGMFs[gmfID][0]] = gmvs[igmv]
        counterGMFs[gmfID] = counterGMFs[gmfID][0] +1

    return uniqueIMT, uniquelocations, setGMFs

def extract_percentile_gmf_xml(gmfFile,percentiles):

    IMTs,exposureLocations,gmfs = parse_nrml_gmf(gmfFile)
    calculate_percentiles_gmfs(percentiles,IMTs,exposureLocations,gmfs)

def save_data_to_csv(data,filename):

    out_file = open(filename,"w")

    for igmv in range(len(data)):
        row = buildRow(data[igmv])
        out_file.write(row[1:]+'\n')
    out_file.close()

def save_mean_std_to_csv(data,filename):

    out_file = open(filename,"w")
    out_file.write('lon,lat,mpga,msa03,msa10,msa30,spga,ssa03,ssa10,ssa30,cpga,csa03,csa10,csa30\n')

    for igmv in range(len(data)):
        singleLoc = data[igmv]
        mPGA, sPGA = calculate_m_std(singleLoc[2],singleLoc[6])
        mSa03, sSa03 = calculate_m_std(singleLoc[3],singleLoc[7])
        mSa10, sSa10 = calculate_m_std(singleLoc[4],singleLoc[8])
        mSa30, sSa30 = calculate_m_std(singleLoc[5],singleLoc[9])
        cPGA = sPGA/mPGA
        cSa03 = sSa03/mSa03
        cSa10 = sSa10/mSa10
        cSa30 = sSa30/mSa30

        out_file.write(str(singleLoc[0])+','+str(singleLoc[1])+','+str(mPGA)+','+str(mSa03)+','+str(mSa10)+','+str(mSa30)+','+str(sPGA)+','+str(sSa03)+','+str(sSa10)+','+str(sSa30)+','+str(cPGA)+','+str(cSa03)+','+str(cSa10)+','+str(cSa30)+'\n')
    out_file.close()

def calculate_m_std(mu,sigma):

    m = math.exp(mu+sigma**2/2)
    std = math.sqrt(math.exp(2*mu+sigma**2)*(math.exp(sigma**2)-1))

    return m, std
