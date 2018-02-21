import math
import numpy as np
from scipy.stats import truncnorm
from scipy import interpolate


def calculate_distance_matrix(exposureLocations):
    noLocations = len(exposureLocations)
    distanceMatrix = np.zeros((noLocations, noLocations))

    for i in range(noLocations):
        for j in range(noLocations):
            distanceMatrix[i, j] = calculate_distance_2points(
                exposureLocations[i], exposureLocations[j])

    return distanceMatrix


def calculate_single_spatial_correlation_matrix(
        distanceMatrix, IMT, correlationType):
    noLocations = len(distanceMatrix)
    correlationMatrix = np.zeros((noLocations, noLocations))
    b = calculate_spatial_length_scale(IMT, 'Vs30clustered')

    for i in range(noLocations):
        for j in range(noLocations):
            if i == j:
                correlationMatrix[i, j] = 1
            else:
                if correlationType == 'no correlation':
                    correlationMatrix[i, j] = 0
                if correlationType == 'full correlation':
                    correlationMatrix[i, j] = 0.99999
                if correlationType == 'spatial':
                    correlationMatrix[i, j] = math.exp(
                        -3 * distanceMatrix[i, j] / b)

    return correlationMatrix


def calculate_spatial_correlation_matrices(
        distanceMatrix, IMTs, correlationType):
    noIMTs = len(IMTs)
    spatialCorrMatrices = []

    for i in range(noIMTs):
        spatialCorrMatrices.append(
            calculate_single_spatial_correlation_matrix(
                distanceMatrix, IMTs[i], correlationType))

    return np.array(spatialCorrMatrices)


def calculate_spatial_covariance_matrices(groundShaking, spatialCorrMatrices):
    # this depends on sPGA, sSa03, sSa10, sSa30
    noIMT = len(spatialCorrMatrices)
    noLocations = len(spatialCorrMatrices[0])
    spatialCovMatrices = []
    for i in range(noIMT):
        tempCovMatrix = np.zeros((noLocations, noLocations))
        for j in range(noLocations):
            for k in range(noLocations):
                tempCovMatrix[j, k] = spatialCorrMatrices[i, j, k] * \
                    groundShaking[j, 2 + noIMT + i] * \
                    groundShaking[k, 2 + noIMT + i]
        spatialCovMatrices.append(tempCovMatrix)

    return np.array(spatialCovMatrices)


def calculate_cross_correlation_matrix(listIMT, correlationType):
    noIMT = len(listIMT)
    crossCorrMatrix = np.zeros((noIMT, noIMT))

    for i in range(noIMT):
        if listIMT[i] == 'PGA':
            T1 = 0.05
        elif listIMT[i][0:2] == 'SA':
            T1 = float(listIMT[i].replace("SA(", "").replace(")", ""))

        for j in range(noIMT):
            if listIMT[j] == 'PGA':
                T2 = 0.05
            elif listIMT[j][0:2] == 'SA':
                T2 = float(listIMT[j].replace("SA(", "").replace(")", ""))

            if i == j:
                crossCorrMatrix[i, j] = 1
            else:

                Tmax = max([T1, T2])
                Tmin = min([T1, T2])

                if Tmin < 0.189:
                    II = 1
                else:
                    II = 0

                if correlationType == 'no correlation':
                    crossCorrMatrix[i, j] = 0
                if correlationType == 'full correlation':
                    crossCorrMatrix[i, j] = 0.99999
                if correlationType == 'cross':
                    crossCorrMatrix[i, j] = 1 - math.cos((math.pi / 2) - (
                        0.359 + 0.163 * II * math.log(Tmin / 0.189)
                    ) * math.log(Tmax / Tmin))

    return crossCorrMatrix


def generate_random_fields_ground_motion(
        IMTs, groundShaking, spatialCovMatrices, crossCorrMatrix,
        siteEffects, noSigmas, noGMFs):
    # groundShaking has shape (N, 11) where
    # 11 = lon lat mPGA mSa03 mSa10 mSa30 sPGA sSa03 sSa10 sSa30 Vs30
    noLocations = spatialCovMatrices.shape[1]
    noIMT = crossCorrMatrix.shape[0]
    L = []
    LLT = []
    Z = []

    for i in range(noIMT):
        L.append(np.linalg.cholesky(spatialCovMatrices[i]))

    L = np.array(L)

    for i in range(noIMT):
        LLTrow = []
        for j in range(noIMT):
            LLTrow.append(
                np.dot(L[i], np.transpose(L[j])) * crossCorrMatrix[i, j])
        for irow in range(len(LLTrow[0])):
            singleLLTrow = np.zeros((int(len(LLTrow) * len(LLTrow[0]))))
            for iL in range(len(LLTrow)):
                singleLLTrow[
                    iL * len(LLTrow[0]):(iL + 1) * len(LLTrow[0])
                ] = LLTrow[iL][irow]
            LLT.append(singleLLTrow)
    LLT = np.array(LLT)

    mu = []
    for i in range(noIMT):
        for j in range(noLocations):
            mu.append(np.ones(noGMFs) * groundShaking[j][i + 2])
    mu = np.array(mu)
    L = (np.linalg.cholesky(LLT))
    Z = truncnorm.rvs(-noSigmas, noSigmas, loc=0, scale=1,
                      size=(noLocations * noIMT, noGMFs))

    gmfs = np.exp(np.dot(L, Z) + mu)

    if siteEffects:  # use vs30 which is the last field
        gmfs = amplifyGMFs(IMTs, groundShaking[:, -1], gmfs) * 0.8

    return gmfs


def amplifyGMFs(IMTs, Vs30s, gmfs):
    noLocations = len(Vs30s)

    for i in range(4):
        if IMTs[i] == 'PGA':
            T = 0.0
        elif IMTs[0:2] == 'SA':
            IMT = IMTs[i]
            T = float(IMT.replace("SA(", "").replace(")", ""))

        for iloc in range(noLocations):
            gmfs[i * noLocations + iloc] = amplify_ground_shaking(
                T, Vs30s[iloc], gmfs[i * noLocations + iloc])

    return gmfs


def calculate_spatial_length_scale(IMT, Vs30Case):
    if IMT == 'PGA':
        T = 0.0
    elif IMT[0:2] == 'SA':
        T = float(IMT.replace("SA(", "").replace(")", ""))

    if T < 1:
        if Vs30Case != 'Vs30clustered':
            b = 8.5 + 17.2 * T
        elif Vs30Case == 'Vs30clustered':
            b = 40.7 - 15.0 * T
    elif T >= 1:
        b = 22.0 + 3.7 * T

    return b


def calculate_distance_2points(point1, point2):
    lon1 = point1[0]
    lon2 = point2[0]
    lat1 = point1[1]
    lat2 = point2[1]

    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * \
        math.cos(lat2) * math.sin(dlon / 2)**2

    distance = 6367 * 2 * math.asin(math.sqrt(a))
    return distance


def amplify_ground_shaking(T, Vs30, IMLs):
    ampFactorsShort = [(760 / Vs30)**0.35,
                       (760 / Vs30)**0.35,
                       (760 / Vs30)**0.25,
                       (760 / Vs30)**0.10,
                       (760 / Vs30)**-0.05,
                       (760 / Vs30)**-0.05]
    ampFactorsMid = [(760 / Vs30)**0.65,
                     (760 / Vs30)**0.65,
                     (760 / Vs30)**0.60,
                     (760 / Vs30)**0.53,
                     (760 / Vs30)**0.45,
                     (760 / Vs30)**0.45]

    if T <= 0.3:
        interpolator = interpolate.interp1d(
            [-1, 0.1, 0.2, 0.3, 0.4, 100],
            ampFactorsShort, kind='linear')
    if T > 0.3:
        interpolator = interpolate.interp1d(
            [-1, 0.1, 0.2, 0.3, 0.4, 100], ampFactorsMid, kind='linear')

    return interpolator(IMLs) * IMLs
