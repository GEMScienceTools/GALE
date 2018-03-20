from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from urllib.request import urlopen
import math
import numpy
import scipy
import zipfile
import urllib
import os


def save_online_data_to_xml(eventID):
    try:
        groundShakingURL = 'http://shakemap.rm.ingv.it/shake/' + \
            str(eventID) + '/download/grid.xml'
        groundShakingData = urlopen(groundShakingURL)
        write_xml(groundShakingData, 'groundShaking' + str(eventID))
    except:
        groundShakingURL = 'http://shakemap.rm.ingv.it/shake/' + \
            str(eventID) + '/download/grid.xml.zip'
        groundShakingData = urlopen(groundShakingURL)
        output = open('groundShaking' + str(eventID) + '.zip', "w")
        output.write(groundShakingData.read())
        output.close()

    try:
        uncertaintyURL = 'http://shakemap.rm.ingv.it/shake/' + \
            str(eventID) + '/download/uncertainty.xml'
        testfile = urllib.URLopener()
        testfile.retrieve(uncertaintyURL, "temp.zip")
        zip_ref = zipfile.ZipFile("temp.zip", 'r')
        zip_ref.extractall('.')
        zip_ref.close()
        os.rename('uncertainty.xml', 'uncertainty' + str(eventID) + '.xml')

    except:
        uncertaintyURL = 'http://shakemap.rm.ingv.it/shake/' + \
            str(eventID) + '/download/uncertainty.xml.zip'
        testfile = urllib.URLopener()
        testfile.retrieve(uncertaintyURL, "temp.zip")
        zip_ref = zipfile.ZipFile("temp.zip", 'r')
        zip_ref.extractall('.')
        zip_ref.close()
        os.rename('uncertainty.xml', 'uncertainty' + str(eventID) + '.xml')


def write_xml(data, filename):
    out_file = open(filename + '.xml', "w")
    for line in data:
        out_file.write(line)
    out_file.close()


def cropShakingData(groundShakingData, typeCrop, limit):
    groundShakingData = numpy.array(groundShakingData)

    lon = groundShakingData[:, 0]
    lat = groundShakingData[:, 1]

    if typeCrop == 'radius':
        epicentre = [scipy.mean(lon), scipy.mean(lat)]

        for iloc in range(len(lon) - 1, 0, -1):
            distance_x = math.fabs(
                epicentre[0] - lon[iloc]) * 111 * math.cos(
                    lat[iloc] * math.pi / 180)
            distance_y = math.fabs(epicentre[1] - lat[iloc]) * 111
            distance = math.sqrt(distance_x**2 + distance_y**2)
            if distance > limit:
                groundShakingData = numpy.delete(groundShakingData, iloc, 0)

    elif typeCrop == 'box':
        for iloc in range(len(lon) - 1, 0, -1):
            if (limit[0] > lon[iloc] or limit[1] < lon[iloc] or
                    limit[2] > lat[iloc] or limit[3] < lat[iloc]):
                groundShakingData = numpy.delete(groundShakingData, iloc, 0)

    print('Total number of locations after cropping')
    print(len(groundShakingData))

    return groundShakingData


def reduceShakingData(groundShakingData, box, res):
    groundShakingData = numpy.array(groundShakingData)
    lon = numpy.linspace(box[0], box[1], (box[1] - box[0]) / res)
    lat = numpy.linspace(box[2], box[3], (box[3] - box[2]) / res)
    reducedGroundShaking = numpy.zeros(
        (int(len(lat) * len(lon)), len(groundShakingData[0])))

    counter = 0
    for ilon in range(len(lon)):
        for ilat in range(len(lat)):
            reducedGroundShaking[counter][0] = lon[ilon]
            reducedGroundShaking[counter][1] = lat[ilat]
            reducedGroundShaking[counter][2:] = findClosestData(
                lon[ilon], lat[ilat], groundShakingData)
            counter = counter + 1

    print('Total number of locations after re-distribution',
          len(reducedGroundShaking))

    return reducedGroundShaking


def findClosestData(lon, lat, groundShakingData):
    dist = numpy.array(
        ((groundShakingData[:, 0] - lon)**2 +
         (groundShakingData[:, 1] - lat)**2)**0.5)
    idx = numpy.where(dist == dist.min())

    return groundShakingData[idx[0], 2:]


def exportLocation(croppedShakingData):
    data = numpy.array(croppedShakingData)
    locations = data[:, 0:2]

    return locations


def parse_xml_data_ID(eventID):
    groundShakingFile = 'groundShaking' + str(eventID) + '.xml'
    uncertaintyFile = 'uncertainty' + str(eventID) + '.xml'
    data, IMTs = parse_xml_data_files(groundShakingFile, uncertaintyFile)

    return data, IMTs


def parse_xml_data_files(groundShakingFile, uncertaintyFile):
    data = []
    file = open(groundShakingFile)
    lines1 = file.readlines()
    file = open(uncertaintyFile)
    lines2 = file.readlines()

    for i in range(len(lines1)):
        line1 = lines1[i].split()
        if isfloat(line1[0]):
            firstFloat1 = i
            break

    for i in range(len(lines2)):
        line2 = lines2[i].split()
        if isfloat(line2[0]):
            firstFloat2 = i
            break

    if len(lines1[firstFloat1].split()) > 9:
        IMTs = ['PGA', 'SA(0.3)', 'SA(1.0)', 'SA(3.0)']
    else:
        IMTs = ['PGA']

    for i in range(len(lines1)):
        line1 = lines1[i].split()
        if isfloat(line1[0]):
            line2 = lines2[i + firstFloat2 - firstFloat1].split()
            lon = float(line1[0])
            lat = float(line1[1])
            if float(line1[2]) > 0:
                sPGA = float(line2[2])
                # sPGA = 0.62
                m = float(line1[2]) / 100
                mPGA = math.log(
                    m**2 / math.sqrt((m**2 * (math.exp(sPGA**2) - 1)) + m**2))
                if len(line1) > 9:
                    sSa03 = float(line2[5])
                    # sSa03 = 0.69
                    m = float(line1[5]) / 100
                    mSa03 = math.log(
                        m**2 / math.sqrt((m**2 * (math.exp(sSa03**2) - 1))
                                         + m**2))

                    sSa10 = float(line2[6])
                    # sSa10 = 0.7
                    m = float(line1[6]) / 100
                    mSa10 = math.log(
                        m**2 / math.sqrt((m**2 * (math.exp(sSa10**2) - 1))
                                         + m**2))

                    sSa30 = float(line2[7])
                    # sSa30 = 0.69
                    m = float(line1[7]) / 100
                    mSa30 = math.log(
                        m**2 / math.sqrt((m**2 * (math.exp(sSa30**2) - 1))
                                         + m**2))

                    Vs30 = float(line1[10])
                    data.append([lon, lat, mPGA, mSa03, mSa10,
                                 mSa30, sPGA, sSa03, sSa10, sSa30, Vs30])
                else:
                    Vs30 = float(line1[7])
                    data.append([lon, lat, mPGA, sPGA, Vs30])
    file.close

    print('Total number of locations with ground shaking:', len(data))
    print('The following IMTs were found:')
    print(IMTs)

    return data, IMTs


def calculateMuSigma(m, sd):

    mu = math.log(m**2 / math.sqrt(sd**2 + m**2))
    sigma = math.sqrt(math.log(sd**2 / m**2 + 1))

    return mu, sigma


def extract_location_exposure(exposureModel):

    assetLocations = []
    uniqueLocations = []

    file = open(exposureModel)
    lines = file.readlines()
    for line in lines:
        if 'location' in line:
            lon = extract_attribute(line, 'lon')
            lat = extract_attribute(line, 'lat')
            assetLocations.append([float(lon), float(lat)])
    file.close

    for i in range(len(assetLocations)):
        if assetLocations[i] not in uniqueLocations:
            uniqueLocations.append(assetLocations[i])

    uniqueLocations = numpy.array(uniqueLocations)

    return uniqueLocations[uniqueLocations[:, 0].argsort()]

#    return uniqueLocations


def extract_IMTs_fragility(fragilityModel):

    IMTs = []
    uniqueIMTs = []
    file = open(fragilityModel)
    lines = file.readlines()
    for line in lines:
        if 'IMT' in line:
            IMT = extract_attribute(line, 'IMT')
            IMTs.append(IMT)
    file.close

    for i in range(len(IMTs)):
        if IMTs[i] not in uniqueIMTs:
            uniqueIMTs.append(IMTs[i])

    return uniqueIMTs


def extract_closest_ground_shaking(exposureLocations, groundShaking):

    groundShaking = numpy.array(groundShaking)
    exposureLocations = numpy.array(exposureLocations)
    closestGroundShaking = numpy.zeros(
        (exposureLocations.shape[0], groundShaking.shape[1]))

    for i in range(exposureLocations.shape[0]):
        location = numpy.array(exposureLocations[i])
        dist_2 = numpy.sum((groundShaking[:, 0:2] - location)**2, axis=1)
        closestGroundShaking[i] = groundShaking[numpy.argmin(dist_2)]
    return closestGroundShaking


def extract_attribute(line, attribute):

    position = line.find(attribute)
    subLine = line[position + len(attribute) + 2:-1]

    return subLine[0:subLine.index('"')]


def isfloat(value):

    try:
        float(value)
        return True
    except ValueError:
        return False


def plotLocations(locations):

    locations = numpy.array(locations)
    box = define_bounding_box(locations[:, 0:2])
    plt.figure(3, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    map = Basemap(llcrnrlon=box["lon_1"], llcrnrlat=box["lat_1"],
                  urcrnrlon=box["lon_2"], urcrnrlat=box["lat_2"],
                  projection='mill', resolution='i')
    x, y = map(locations[:, 0], locations[:, 1])

    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.drawstates(linewidth=0.25)
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='white', lake_color='aqua')
    plt.scatter(x, y, 5, marker='o', color='k', zorder=4)

    if box["lat_length"] + 2 < 4:
        parallels = numpy.arange(0., 81, 0.5)
    else:
        parallels = numpy.arange(0., 81, 1.0)
    map.drawparallels(parallels, labels=[True, False, True, False])
    if box["lon_length"] + 2 < 4:
        meridians = numpy.arange(0., 360, 0.5)
    else:
        meridians = numpy.arange(0., 360, 1.0)
    map.drawmeridians(meridians, labels=[True, False, False, True])

    plt.show()


def define_bounding_box(locations):

    box = {"lon_0": None, "lat_0": None, "lon_1": None,
           "lat_1": None, "lon_2": None, "lat_2": None}

    maxCoordinates = locations.max(axis=0)
    minCoordinates = locations.min(axis=0)
    lengthLon = abs(maxCoordinates[0] - minCoordinates[0])
    lengthLat = abs(maxCoordinates[1] - minCoordinates[1])

    box["lon_0"] = (minCoordinates[0] + maxCoordinates[0]) / 2
    box["lat_0"] = (minCoordinates[1] + maxCoordinates[1]) / 2
    box["lon_1"] = minCoordinates[0] - 1
    box["lat_1"] = minCoordinates[1] - 1
    box["lon_2"] = maxCoordinates[0] + 1
    box["lat_2"] = maxCoordinates[1] + 1
    box["lon_length"] = lengthLon
    box["lat_length"] = lengthLat

    return box
