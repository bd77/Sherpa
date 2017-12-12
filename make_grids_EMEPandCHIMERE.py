'''
Created on Jun 28, 2017

creates the CHimere and EMEP grids in the database

CREATE TABLE sherpa.grids
(
  id serial NOT NULL,
  grid_name character varying(25),
  cell_id integer,
  the_cell geometry,
  CONSTRAINT pkey_id PRIMARY KEY (id)
)

@author: degraba
'''

import psycopg2

host_name = "localhost"
database_name = "sherpa"
user_name = "postgres"
pw = "root"
schema_name = 'sherpa'
table_name = 'grids'


# database connection
try:
    conn = psycopg2.connect(host=host_name, database=database_name, user=user_name, password=pw);
except:
    print("I am unable to connect to the database")
cur = conn.cursor()

# remove the old grid from the table
if 1 == 1:
    cur.execute('DELETE FROM sherpa.grids;')
    conn.commit()    

# grid properties
grid_name_lst = ['chimere', 'emep']
xmin_dict = {'chimere': -10.5, 'emep': -30.0}  
xmax_dict = {'chimere': 37.5, 'emep': 90.0}    
ymin_dict = {'chimere': 34.0, 'emep': 30.0} 
ymax_dict = {'chimere': 62, 'emep': 82.0}
dlon_dict = {'chimere': 0.125, 'emep': 0.1}
dlat_dict = {'chimere': 0.0625, 'emep': 0.1}
srid_WGS84 = 4326

# loop over all grids (Chimere and EMEP)
for grid_name in grid_name_lst:
    xmin = xmin_dict[grid_name]
    xmax = xmax_dict[grid_name]
    ymin = ymin_dict[grid_name]
    ymax = ymax_dict[grid_name]
    dlon = dlon_dict[grid_name]                # delta lon in degrees
    dlat = dlat_dict[grid_name]                # delta lat in degrees

    # create all cells of the grid from south to north and west to east
    # initialize cell id and coordinates
    cell_id = 1
    x_west = xmin
    x_east = xmin + dlon
    y_south = ymin
    y_north = ymin + dlat
    
    while x_east <= xmax:
        
        while y_north <= ymax:
            
            # wkt for the polygon of the grid cell
            polygon_wkt = "ST_GeomFromText('POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))', %d)" \
            % (x_west, y_south, x_east, y_south, x_east, y_north, x_west, y_north, x_west, y_south, srid_WGS84)
            query = "INSERT INTO %s.%s(grid_name, cell_id, the_cell) VALUES ('%s', %d, %s);" % (schema_name, table_name, grid_name, cell_id, polygon_wkt)
            # print(query)
            cur.execute(query)
            conn.commit()
            cell_id += 1
            # go one cell north
            y_south += dlat
            y_north += dlat
    
        # go one cell east an go back south
        x_east += dlon
        x_west += dlon
        y_south = ymin
        y_north = ymin + dlat
    
    # create table ifdm_ia.emission_gridding as select ST_Intersection(the_road, the_cell) from ifdm.lijnbroninvoer_antwerp, (select * from ifdm_ia.grids where cell_id=119) g