import csv

def dms_to_decimal(dms_str):
    parts = dms_str.split()
    deg, min, sec = map(int, parts)

    # Check if degrees are negative
    if deg < 0:
        return deg - min / 60 - sec / 3600
    else:
        return deg + min / 60 + sec / 3600

def remove_duplicates_get_indices(lst):
    seen = set()
    seen_add = seen.add
    no_duplicates_list = []
    duplicates_indices = []
    for i, item in enumerate(lst):
        if item not in seen:
            no_duplicates_list.append(item)
            seen_add(item)
        else:
            duplicates_indices.append(i)
    return no_duplicates_list, duplicates_indices

def remove_elements_from_indices(lst, indices):
    return [i for j, i in enumerate(lst) if j not in indices]

def convert_to_floats(mixed_list, default=None):
    float_list = []
    for i in mixed_list:
        try:
            float_list.append(float(str(i).replace(" ", "")))
        except ValueError:
            float_list.append(default)
    return float_list

stations=[]
lats=[]
lons=[]
ids=[]

# Read this file first since it includes actual IDs and we want to keep those in
# the event of dupes
with open('all_dcp_defs.txt', 'r') as f:
    reader = csv.reader(f, delimiter='|')
    for row in reader:
        decimal = dms_to_decimal(row[5])
        stations.append(row[1])
        lats.append(dms_to_decimal(row[5]))
        lons.append(dms_to_decimal(row[6]))
        ids.append(row[0])

# Open the file
with open('NWSLI20161104.TXT', 'r') as file:
    # Create a CSV reader, specify the delimiter as "|"
    reader = csv.reader(file, delimiter='|')

    # Process each row in the file
    for row in reader:
        lats.append(row[-3])
        lons.append(row[-2])
        stations.append(row[1])
        ids.append('XXX')
##



lats = convert_to_floats(lats)
lons = convert_to_floats(lons)

stations_no_duplicates, duplicates_indices = remove_duplicates_get_indices(stations)

lats_no_duplicates = remove_elements_from_indices(lats, duplicates_indices)
lon_no_duplicates = remove_elements_from_indices(lons, duplicates_indices)
ids_no_duplicates = remove_elements_from_indices(ids, duplicates_indices)

# ensure all lists are of the same length
assert len(stations_no_duplicates) == len(lats_no_duplicates) == len(lon_no_duplicates) == len(ids_no_duplicates)

# open the file in write mode
with open('processed-stations.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    # write the headers
    # CMZ, no header needed since passing into NCL
    #writer.writerow(['Station ID', 'ID', 'Latitude', 'Longitude'])

    # write the data
    for station_id, idname, latitude, longitude in zip(stations_no_duplicates, ids_no_duplicates, lats_no_duplicates, lon_no_duplicates):
        writer.writerow([station_id, idname, latitude, longitude])

