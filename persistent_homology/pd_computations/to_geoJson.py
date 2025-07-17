import csv
import json

def csv_to_geojson(csv_file, output_file=None):
    geojson = {"type": "FeatureCollection", "features": []}
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        for i, row in enumerate(reader, 2):
            try:
                if not row.get('lat') or not row.get('lng'):
                    print(f"Row {i}: Missing lat or lng")
                    continue
                    
                lat = float(row['lat'])
                lng = float(row['lng'])
                
                if not (-90 <= lat <= 90) or not (-180 <= lng <= 180):
                    print(f"Row {i}: Invalid coordinates")
                    continue
                
                if not row.get('name') or not row.get('address'):
                    print(f"Row {i}: Missing name or address")
                    continue
                
                feature = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": [lng, lat]
                    },
                    "properties": {
                        "name": row['name'].strip(),
                        "address": row['address'].strip()
                    }
                }
                geojson["features"].append(feature)
                
            except (ValueError, KeyError) as e:
                print(f"Row {i}: Error - {e}")
                continue
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(geojson, f, indent=2)
    
    print(f"Processed {len(geojson['features'])} locations")
    return geojson

# Usage
data_all = csv_to_geojson('persistent_homology/coordinates_all.csv', 'output_all.geojson')
data_fqhc = csv_to_geojson('persistent_homology/coordinates_fqhc.csv', 'output_fqhc.geojson')