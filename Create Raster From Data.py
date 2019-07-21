# Import math
import math

# Import arcpy module
import arcpy

# Environment Specs
arcpy.env.overwriteOutput = True

# Interpolation Function
def Barnes_Interp(x_short, y_short, stats_at_date, new_x, new_y):

    gamma = 0.2

    del_n_c = 0
    h = len(x_short)
    for m in range(0, h):
        d = []
        for n in range(0, h):
            d.append(math.sqrt((x_short[n]-x_short[m])**2 + (y_short[n]-y_short[m])**2))

        d.sort()

        del_n_c += d[1]/h

    area = (max(x_short)-min(x_short))*(max(y_short)-min(y_short))

    del_n_r = math.sqrt(area)*((1+math.sqrt(h))/(h-1))

    del_n = ((del_n_r - del_n_c)/2) + del_n_c

    k = 5.052*(2*del_n/math.pi)**2

    g0 = []
    new_stats = []

    for m in range(0, len(new_y)):
        d_0 = []
        weight = []
        num = 0
        den = 0
        val_sum = 0

        for n in range(0, h):

            d_0.append(math.sqrt((new_x[m] - x_short[n])**2 + (new_y[m] - y_short[n])**2))
            weight.append(math.exp(-1*(d_0[n]**2)/k))

            if stats_at_date[n] is None:
                continue

            if isinstance(stats_at_date[n], basestring):
                stats_at_date[n] = float(stats_at_date[n])
            if isinstance(weight[n], basestring):
                weight[n] = float(weight[n])

            num = num + stats_at_date[n]*weight[n]
            den = den + weight[n]

        g0.append(num/den)

        for n in range(0, h):
            if stats_at_date[n] is None:
                continue
            val_sum = val_sum + ((stats_at_date[n] - g0[m])*(math.exp(-1*(d_0[n]**2)/(gamma*k))))

        g1 = g0[m] + val_sum
        if g1 < 0:
            g1 = 0
        new_stats.append(g1)

    return new_stats


# Script Arguments
Data_Table = arcpy.GetParameterAsText(0)

Param_Name = arcpy.GetParameterAsText(1)

Location_Data = arcpy.GetParameterAsText(2)

Output_Join_Field = arcpy.GetParameterAsText(3)
if Output_Join_Field == '#' or not Output_Join_Field:
    Output_Join_Field = "Station"  # provide a default value if unspecified

Join_Fields = arcpy.GetParameterAsText(4)

Interp_Table = arcpy.GetParameterAsText(5)

Del_Frac = arcpy.GetParameterAsText(6)
if Del_Frac == '#' or not Del_Frac:
    Del_Frac = 20

Del_Frac = float(Del_Frac)

Cellsize = arcpy.GetParameterAsText(7)

Cell_assignment_type = arcpy.GetParameterAsText(8)
if Cell_assignment_type == '#' or not Cell_assignment_type:
    Cell_assignment_type = "MEAN" # provide a default value if unspecified

RASTER_Name = arcpy.GetParameterAsText(9)

X_header = arcpy.GetParameterAsText(10)
Y_header = arcpy.GetParameterAsText(11)

# Spatial Reference
Spatial_Reference = arcpy.Describe(Location_Data).spatialReference
Spat_Ref = Spatial_Reference.exportToString()

# Reading data from table
arcpy.AddMessage("Reading Data from Table")
header_list = []
fields = arcpy.ListFields(Data_Table)

for x in fields:
    header_list.append(str(x.name))

date_list = []
stat_list = [[] for Null in range(0, len(header_list)-2)]
name_list = [[] for Null in range(0, len(header_list)-2)]

with arcpy.da.SearchCursor(Data_Table, "*") as cur:
    for row in cur:

        date_list.append(row[1])

        for x in range(2, len(header_list)):
            stat_list[x-2].append(row[x])

            name_list[x-2].append(header_list[x])

        # End for
del cur

# Creating New Table
arcpy.AddMessage("Creating New Table")
idx = Data_Table.rfind("\\")
Table_name = "Test_Table_1"
pre_Work_Space = Data_Table[:idx]

idx2 = pre_Work_Space.rfind("\\")
Work_Space = pre_Work_Space[:idx2]

arcpy.env.workspace = Work_Space

New_Table = pre_Work_Space + "\\" + Table_name

arcpy.CreateTable_management(pre_Work_Space, Table_name)
arcpy.AddField_management(New_Table, "Date", "DATE")
arcpy.AddField_management(New_Table, Param_Name, "FLOAT")
arcpy.AddField_management(New_Table, "Location", "TEXT")

# Adding Data to Table
arcpy.AddMessage("Adding Data to New Table")
rows = arcpy.InsertCursor(New_Table)

for x in range(2, len(header_list)):

    for y in range(0, len(date_list)):
        row = rows.newRow()
        row.setValue("Date", date_list[y])
        row.setValue(Param_Name, stat_list[x-2][y])
        row.setValue("Location", name_list[x-2][y])

        rows.insertRow(row)

del rows

# Joining Location Data to Table
arcpy.AddMessage("Joining Table Data to Location Data")
arcpy.JoinField_management(New_Table, "Location", Location_Data, Output_Join_Field, Join_Fields)

# Reading Location Data From Table
arcpy.AddMessage("Getting Location Data From Table")

fields = arcpy.ListFields(New_Table)

X_header_test = str(fields[4].name)

if X_header == X_header_test:
    col_x = 4
    col_y = 5
else:
    col_x = 5
    col_y = 4

X_list = [[] for Null in range(0, len(header_list)-2)]
Y_list = [[] for Null in range(0, len(header_list)-2)]

x = 0
y = 1

with arcpy.da.SearchCursor(New_Table, "*") as cur:
    for row in cur:
        if y > len(date_list):
            x += 1
            y = 1

        X_list[x].append(row[col_x])
        Y_list[x].append(row[col_y])

        y += 1

# Determining New Grid Points
arcpy.AddMessage("Determining New Grid Point Locations")
X_short_list = []
Y_short_list = []

delta_X = []
delta_Y = []
for x in range(0, len(header_list)-2):
    X_short_list.append(X_list[x][0])
    Y_short_list.append(Y_list[x][0])

    if x > 1:
        del_x = abs(X_list[x][0] - X_list[x-1][0])
        del_y = abs(Y_list[x][0] - Y_list[x - 1][0])
        if del_x > 0:
            delta_X.append(del_x)
        if del_y > 0:
            delta_Y.append(del_y)

min_X = min(X_short_list)
max_X = max(X_short_list)
min_delta_X = min(delta_X)/Del_Frac

min_Y = min(Y_short_list)
max_Y = max(Y_short_list)
min_delta_Y = min(delta_Y)/Del_Frac

New_X = [min_X]
New_Y = [min_Y]

x = 0
while New_X[x] < max_X:
    New_X.append(New_X[x] + min_delta_X)
    x = x + 1

x = 0
while New_Y[x] < max_Y:
    New_Y.append(New_Y[x] + min_delta_Y)
    x = x + 1

New_X_long = []
New_Y_long = []

for x in range(0, len(New_X)):

    for y in range(0, len(New_Y)):
        New_X_long.append(New_X[x])
        New_Y_long.append(New_Y[y])

# Interpolating New Grid Values
arcpy.AddMessage("Interpolating New Grid Values")

New_stat_list = [[] for Null in range(0, len(New_X_long))]

for x in range(0, len(date_list)):

    stat_at_date = []

    for y in range(0, len(X_short_list)):
        stat_at_date.append(stat_list[y][x])


    New_stat_at_date = Barnes_Interp(X_short_list, Y_short_list, stat_at_date, New_X_long, New_Y_long)

    for y in range(0, len(New_X_long)):
            New_stat_list[y].append(New_stat_at_date[y])

# Creating New Table of Interpolated Values
arcpy.AddMessage("Creating Table of Interpolated Values")

idx3 = Interp_Table.rfind("\\")
Interp_Table_Name = Interp_Table[idx3+1:]

arcpy.CreateTable_management(pre_Work_Space, Interp_Table_Name)
arcpy.AddField_management(Interp_Table, "Date", "DATE")
arcpy.AddField_management(Interp_Table, Param_Name, "DOUBLE")
arcpy.AddField_management(Interp_Table, X_header, "DOUBLE")
arcpy.AddField_management(Interp_Table, Y_header, "DOUBLE")

# Adding Data to Interp Table
arcpy.AddMessage("Adding Data to Table of Interpolated Values")

rows = arcpy.InsertCursor(Interp_Table)

for x in range(0, len(New_X_long)):

    for y in range(0, len(date_list)):
        row = rows.newRow()
        row.setValue("Date", date_list[y])
        row.setValue(Param_Name, New_stat_list[x][y])
        row.setValue(X_header, New_X_long[x])
        row.setValue(Y_header, New_Y_long[x])

        rows.insertRow(row)

del rows

# Deleting Test Table
arcpy.AddMessage("Deleting Joined Table")
arcpy.Delete_management(New_Table)

# Process: Make XY Event Layer
arcpy.AddMessage("Creating Point Feature Class From Table")

Layer_Name = Interp_Table_Name + "_Points"
Layer_Path = pre_Work_Space + "\\" + Layer_Name

arcpy.MakeXYEventLayer_management(Interp_Table, X_header, Y_header, Layer_Path, Spat_Ref, Param_Name)
arcpy.FeatureClassToFeatureClass_conversion(Layer_Path, pre_Work_Space, Layer_Name)

# Process: Point to Raster
arcpy.AddMessage("Creating the Raster")
arcpy.PointToRaster_conversion(Layer_Path, Param_Name, RASTER_Name, Cell_assignment_type, "NONE", Cellsize)

# Deleting Stuff
arcpy.AddMessage("Deleting Unnecessary Feature Classes")
arcpy.Delete_management(Interp_Table)
arcpy.Delete_management(Layer_Path)
