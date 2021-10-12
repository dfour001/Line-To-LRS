import arcpy, pandas as pd

inputEventLayer = r'C:\Users\daniel.fourquet\Desktop\Tasks\XD Speed Limits for JJ\data.gdb\requiredXDsEventTable'
idField = 'id'
rte_nmField = 'rte_nm'
begin_mpField = 'begin_mp'
end_mpField = 'end_mp'

outputEventCSV = r'C:\Users\daniel.fourquet\Documents\ArcGIS\Projects\Scrap\outputFlipped.csv'

overlapLRS = r'C:\Users\daniel.fourquet\Documents\ArcGIS\LRS.gdb\SDE_VDOT_RTE_OVERLAP_LRS'


outputEvents = []

def add_to_event_table(id, rte_nm, begin_mp, end_mp):
    event = {
        'id': id,
        'rte_nm': rte_nm,
        'begin_mp': begin_mp,
        'end_mp': end_mp
    }

    outputEvents.append(event)


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
    beginMP = get_mp_from_point(RouteGeom, beginPt)

    endPt = testGeom.lastPoint
    endMP = get_mp_from_point(RouteGeom, endPt)

    return round(beginMP, 3), round(endMP, 3)


# Create a dictionary of opposite direction routes
inputRoutes = set([row[0] for row in arcpy.da.SearchCursor(inputEventLayer, 'rte_nm')])
oppRteDict = {}
with arcpy.da.SearchCursor(overlapLRS, ['RTE_NM', 'RTE_OPPOSITE_DIRECTION_RTE_NM']) as cur:
    for rte_nm, opp_rte_nm in cur:
        if rte_nm in inputRoutes:
            oppRteDict[rte_nm] = opp_rte_nm


# For each record in input layer, if the begin_mp > end_mp, move to the opposite route.  Otherwise, keep the same
with arcpy.da.SearchCursor(inputEventLayer, [idField, rte_nmField, begin_mpField, end_mpField, 'SHAPE@']) as cur:
    for id, rte_nm, begin_mp, end_mp, geom in cur:
        try:
            if begin_mp < end_mp or begin_mp == end_mp:
                add_to_event_table(id, rte_nm, begin_mp, end_mp)
            else:
                new_rte_nm = oppRteDict[rte_nm]
                begin_mp, end_mp = get_msr(geom, overlapLRS, new_rte_nm)
                add_to_event_table(id, new_rte_nm, begin_mp, end_mp)
        except:
            print(f'Error processing {id}')
            add_to_event_table(id, rte_nm, begin_mp, end_mp)

df = pd.DataFrame(outputEvents)
df.to_csv(outputEventCSV, index=False)