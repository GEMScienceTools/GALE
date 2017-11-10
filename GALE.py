import urllib2
import math
import numpy
import scipy
import inputProcessing
import groundShakingCalculation
import outputProcessing
import time

#2957481

eventID = 1

print 'ACCESSING INGV SHAKEMAP DATABASE...'
#inputProcessing.save_online_data_to_xml(eventID)

print 'SAVING INGV DATA INTO XML...'
groundShakingData, IMTs = inputProcessing.parse_INGV_xml_data(eventID)
#IMTs = ['PGA', 'SA(0.3)']

outputProcessing.save_data_to_csv(groundShakingData,'rawtnepal.csv')

#croppedShakingData = inputProcessing.cropShakingData(groundShakingData,'box',[10.5,12,44.5,45.5])

#reducedShakingData = inputProcessing.reduceShakingData(groundShakingData,[10.5,12,44.2,45.8],0.02)

#outputProcessing.save_mean_std_to_csv(reducedShakingData,'statdata.csv')

#locations =  inputProcessing.exportLocation(reducedShakingData)
locations = inputProcessing.extract_location_exposure('exposure_test.xml')
print locations

print 'CALCULATING CLOSEST GROUND SHAKING...'
closestGroundShaking = inputProcessing.extract_closest_ground_shaking(locations,groundShakingData)
outputProcessing.save_data_to_csv(closestGroundShaking,'closestnepal.csv')
print closestGroundShaking

print 'CALCULATING DISTANCE BETWEEN LOCATIONS...'
distanceMatrix = groundShakingCalculation.calculate_distance_matrix(locations)
print distanceMatrix

print 'CALCULATING SPATIAL CORRELATION MATRIX...'
spatialCorrMatrices = groundShakingCalculation.calculate_spatial_correlation_matrices(distanceMatrix,IMTs,'spatial')
print spatialCorrMatrices

print 'CALCULATING SPATIAL COVARIANCE MATRIX...'
spatialCovMatrices = groundShakingCalculation.calculate_spatial_covariance_matrices(closestGroundShaking,spatialCorrMatrices)
print spatialCovMatrices

print 'CALCULATING CROSS-CORRELATION MATRIX...'
crossCorrMatrix = groundShakingCalculation.calculate_cross_correlation_matrix(IMTs,'cross')
print crossCorrMatrix

#print 'GENERATING SPATIALLY-CROSS-CORRELATED RANDOM FIEDS OF GROUND MOTION...'
gmfs = groundShakingCalculation.generate_random_fields_ground_motion(closestGroundShaking,spatialCovMatrices,crossCorrMatrix,2000)
#print 'SAVING GROUND MOTION FIELDS IN NRML FORMAT..'

outputProcessing.save_gmfs_to_csv(eventID,locations,gmfs)
#print IMTs
#print exposureLocations
#print gmfs
outputProcessing.save_gmfs_to_xml(eventID,IMTs,locations,gmfs)
#outputProcessing.calculate_percentiles_gmfs([16,50,84],IMTs,exposureLocations,gmfs)