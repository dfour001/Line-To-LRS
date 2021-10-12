import arcpy
import numpy
import json
from scipy.spatial.distance import directed_hausdorff
import pandas as pd
from datetime import datetime
import logging

# Search Distance in Meters
d = 50

# LRS Layer
lrs = r'C:\Users\daniel.fourquet\Documents\ArcGIS\LRS.gdb\SDE_VDOT_RTE_MASTER_LRS'

# LRS Feature Layer - prime direction and no ramps:
# arcpy.MakeFeatureLayer_management(lrs, 'memory\lrs_noRamps', "RTE_DIRECTION_CD = RTE_PRIME_DIRECTION_CD and rte_nm not like '%RMP%'")
# lrs_noRamps = 'memory\lrs_noRamps'

# No Ramps:
# arcpy.MakeFeatureLayer_management(lrs, 'memory\lrs_noRamps', "rte_nm not like '%RMP%'")
# lrs_noRamps = 'memory\lrs_noRamps'

# All LRS
arcpy.MakeFeatureLayer_management(lrs, 'lrsAll', "RTE_DIRECTION_CD = RTE_PRIME_DIRECTION_CD")
lrs_noRamps = "lrsAll"


# Only Ramps:
arcpy.MakeFeatureLayer_management(lrs, 'memory\lrs_onlyRamps', "rte_nm like '%RMP%'")
lrs_onlyRamps = 'memory\lrs_onlyRamps'

# LRS Intersections
Intersections = r'C:\Users\daniel.fourquet\Documents\ArcGIS\LRS.gdb\SDE_VDOT_INTERSECTION_W_XY'
lyrIntersections = "lyrIntersections"
arcpy.MakeFeatureLayer_management(Intersections, lyrIntersections)

# Input lines file path
inputLines = r'C:\Users\daniel.fourquet\Desktop\Tasks\XD Speed Limits for JJ\data.gdb\RequiredXDs'

# Output CSV file path
outputPath = r'C:\Users\daniel.fourquet\Documents\GitHub\Line-To-LRS\testData\SpeedLimitTest.csv'



# Output events
events = []

# Set up logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) # Set the debug level herec
fileHandler = logging.FileHandler(r'C:\Users\daniel.fourquet\Documents\GitHub\Line-To-LRS\testData\test.log', mode='w')
log.addHandler(fileHandler)

arcpy.env.overwriteOutput = True


# Spatial reference
sr = arcpy.SpatialReference(3857)


def get_array_from_geom(geom):
    esriJSON = json.loads(geom.JSON)
    geomPoints = []
    for path in esriJSON['paths']:
        for point in path:
            xy = [point[0], point[1]]
            geomPoints.append(xy)
    
    array = numpy.array(geomPoints)
    return array


def get_angle_distance(geom):
    firstPoint = arcpy.PointGeometry(geom.firstPoint, sr)
    lastPoint = arcpy.PointGeometry(geom.lastPoint, sr)
    
    angle, distance = firstPoint.angleAndDistanceTo(lastPoint)
    
    if angle < 0:
        angle = angle + 180
        
    return angle, distance


def clip_lrs(geom, lrs):
    InputBuffer = geom.buffer(d)
    lrsClip = arcpy.analysis.Clip(lrs, InputBuffer, "memory/lrsClip") 
    return lrsClip


def get_hausdorff(testGeom, routeGeom):
    testArray = get_array_from_geom(testGeom)
    routeArray = get_array_from_geom(routeGeom)
    hausdorff = directed_hausdorff(testArray, routeArray)[0]
    
    return hausdorff


def move_to_closest_int(testGeom, testDistance=30):
    testPointGeom = arcpy.PointGeometry(testGeom, arcpy.SpatialReference(3857))
    arcpy.SelectLayerByLocation_management(lyrIntersections, "INTERSECT", testPointGeom, testDistance)
    intersections = [row[0] for row in arcpy.da.SearchCursor(lyrIntersections, "SHAPE@")]
    if len(intersections) == 0:
        log.debug(f"No intersections within {testDistance}m distance.  Returning testGeom.")
        return testGeom
    
    if len(intersections) == 1:
        log.debug(f"One intersection within {testDistance}m distance.  Returning intersection at {intersections[0].firstPoint.X}, {intersections[0].firstPoint.Y}.")
        return intersections[0].firstPoint

    log.debug(f"{len(intersections)} intersections found.  Returning closest intersection.")
    closestInt = intersections[0]
    closestIntDist = testDistance
    for intersection in intersections:
        dist = testPointGeom.distanceTo(intersection)
        if dist < closestIntDist:
            closestInt = intersection
            closestIntDist = dist
    
    log.debug(f"Returning intersection at {closestInt.firstPoint.X}, {closestInt.firstPoint.Y}")
    return closestInt.firstPoint

    



def get_msr(testGeom, lrs, rte_nm):
    with arcpy.da.SearchCursor(lrs, "SHAPE@", "RTE_NM = '{}'".format(rte_nm)) as cur:
        for row in cur:
            RouteGeom = row[0]

    def get_mp_from_point(route, point):
        rteMeasure = route.measureOnLine(point)
        rtePosition = route.positionAlongLine(rteMeasure)
        mp = rtePosition.firstPoint.M
        return mp

    beginPt = testGeom.firstPoint
    beginClosestPt = move_to_closest_int(beginPt)
    

    endPt = testGeom.lastPoint
    endClosestPt = move_to_closest_int(endPt)

    if beginClosestPt.X != endClosestPt.X and beginClosestPt.Y != endClosestPt.Y:
        log.debug("Begin/End Closest Point are not equal")
        beginPt = beginClosestPt
        endPt = endClosestPt
    else:
        log.debug("BEGIN/END CLOSEST POINT ARE EQUAL")

    beginMP = get_mp_from_point(RouteGeom, beginPt)
    endMP = get_mp_from_point(RouteGeom, endPt)

    return round(beginMP, 3), round(endMP, 3)


def normalize_scores(scores):
    normalized = []
    
    minScore = min(scores)
    maxScore = max(scores)
    for score in scores:
        normScore = round((score - minScore) / (maxScore - minScore), 2) * 100
        normalized.append(normScore)
    
    return normalized


def normalize_score(score, minScore, maxScore):
    normScore = round((score - minScore) / (maxScore - minScore), 2) * 100
    
    return normScore


def add_event(testGeom, ID, rte_nm):
    beginMP = None
    endMP = None
    if rte_nm:
        beginMP, endMP = get_msr(testGeom, lrs, rte_nm)
        
    outputEvent = {
        'id': ID,
        'rte_nm': rte_nm,
        'begin_mp': beginMP,
        'end_mp': endMP
    }
            
    events.append(outputEvent)


class PotentialRouteMatch:
    def __init__(self, rte_nm, hausdorff, angleDiff, distanceDiff, score_rte_nbr, matchScore):
        self.rte_nm = rte_nm
        self.hausdorff = hausdorff
        self.angleDiff = angleDiff
        self.score_rte_nbr = score_rte_nbr
        self.matchScore = round(matchScore, 2)
        self.matchScoreNormalized = None
        
    def normalize(self, allScores):
        self.matchScoreNormalized = normalize_score(self.matchScore, min(allScores), max(allScores))
        log.debug(f'{self.rte_nm}: {self.matchScore} => {self.matchScoreNormalized}')
        
    def __repr__(self):
        return f'\n{self.rte_nm}:\n  hausdorff: {self.hausdorff}\n  andleDiff: {self.angleDiff}\n  score_rte_nbr: {self.score_rte_nbr}\n\n  MATCH SCORE: {self.matchScore}\n\n'
    


def create_events(testGeom, ID, lrs, rte_nbr=-1, ):
    nowProcessing = f'== {ID} =='
    log.debug('='*len(nowProcessing))
    log.debug(nowProcessing)
    log.debug('='*len(nowProcessing))
    testGeom = testGeom.projectAs(sr)
    lrsClip = clip_lrs(testGeom, lrs)
    
    try:
        rte_nbr = int(rte_nbr)
    except:
        rte_nbr = -1

    potentialRouteMatches = []

    with arcpy.da.SearchCursor(lrsClip, ['RTE_NM', 'SHAPE@', 'RTE_NBR']) as cur:
        for row in cur:
            rte_nm = row[0]
            routeGeom = row[1]
            lrs_rte_nbr = row[2]

            hausdorff1 = get_hausdorff(testGeom, routeGeom)
            hausdorff2 = get_hausdorff(routeGeom, testGeom)
            # hausdorff = min(hausdorff1, hausdorff2)
            hausdorff = (hausdorff1 + hausdorff2) / 2

            
            testGeomAngle, testGeomDistance = get_angle_distance(testGeom)
            routeGeomAngle, routeGeomDistance = get_angle_distance(routeGeom)
            
            angleDiff = abs(testGeomAngle - routeGeomAngle)
            if angleDiff > 50:
                angleDiff = 200
                
            distanceDiff = abs(testGeomDistance - routeGeomDistance)
            
            
            if lrs_rte_nbr and rte_nbr == lrs_rte_nbr:
                score_rte_nbr = 0.75
            else:
                score_rte_nbr = 1
                
                    
            
            matchScore = (hausdorff + angleDiff) * score_rte_nbr
            
            potentialRouteMatch = PotentialRouteMatch(rte_nm, hausdorff, angleDiff, distanceDiff, score_rte_nbr, matchScore)
            potentialRouteMatches.append(potentialRouteMatch)


    potentialRouteMatches = sorted(potentialRouteMatches, key=lambda x: x.matchScore)
    log.debug(potentialRouteMatches)
    
    if len(potentialRouteMatches) > 1:
        scores = [match.matchScore for match in potentialRouteMatches]
        log.debug('\nScores:')
        log.debug(scores)
    
        for potentialRouteMatch in potentialRouteMatches:
            potentialRouteMatch.normalize(scores)
    
        scoresNormalized = [match.matchScoreNormalized for match in potentialRouteMatches]
        log.debug('\nScores Normalized:')
        log.debug(scoresNormalized)
    elif len(potentialRouteMatches) == 0:
        add_event(testGeom, ID, None)
    else:
        potentialRouteMatches[0].matchScoreNormalized = 1

    ramp = False
    if 'RMP' in potentialRouteMatches[0].rte_nm:
        # This event most likely represents a ramp, so only ramp routes should be matched
        ramp = True
    
    for potentialRouteMatch in potentialRouteMatches:
        if potentialRouteMatch.matchScoreNormalized <= 10:
            if ramp:
                if 'RMP' in potentialRouteMatch.rte_nm:
                    add_event(testGeom, ID, potentialRouteMatch.rte_nm)
            else:
                add_event(testGeom, ID, potentialRouteMatch.rte_nm)
            



if __name__ == '__main__':

    # Configuration for XD-LRS conflation:
    startTime = datetime.now()
    print(f'Start time: {startTime}')
    XDnonRamps = [row[0] for row in arcpy.da.SearchCursor(inputLines, ["XDSegID", "SlipRoad", "SHAPE@LENGTH"]) if row[1] != '1' or (row[1] == '1' and row[2] > 0.009)]
    XDramps = [row[0] for row in arcpy.da.SearchCursor(inputLines, ["XDSegID", "SlipRoad", "SHAPE@LENGTH"]) if row[1] != '0' and row[2] < 0.009]

    # Process Non-ramps
    print('\n=== Processing Non-Ramps ===')
    log.debug('\n=== Processing Non-Ramps ===')
    for i, testID in enumerate(XDnonRamps):
        try:
            testGeom, RoadNumber = [(row[0], row[1]) for row in arcpy.da.SearchCursor(inputLines, ["SHAPE@", 'RoadNumber'], f"XDSegID = {testID}")][0]
            create_events(testGeom, testID,  lrs_noRamps, RoadNumber)
        except Exception as e:
            print(e)
            print(f'Error processing {testID}')

        percent = f'{i / len(XDnonRamps) * 100}%'
        print(round(percent))


    # Process Ramps
    print('\n=== Processing Ramps ===')
    log.debug('\n=== Processing Ramps ===')
    for i, testID in enumerate(XDramps):
        try:
            testGeom, RoadNumber = [(row[0], row[1]) for row in arcpy.da.SearchCursor(inputLines, ["SHAPE@", 'RoadNumber'], f"XDSegID = {testID}")][0]
            create_events(testGeom, testID, lrs_onlyRamps, RoadNumber)
        except:
            print(f'Error processing {testID}')

        percent = f'{i / len(XDramps) * 100}%'
        print(round(percent))

    outputDF = pd.DataFrame(events)
    outputDF.to_csv(outputPath, index=False)

    print(f'\n\nComplete: {datetime.now() - startTime}')