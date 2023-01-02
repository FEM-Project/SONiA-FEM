# Quarter of Plate 2 Elements:
function quarterPlate_2el()
    # Connectivity Matrix
    conn = [1 4 1 1 2 3 4 
        2 4 1 4 3 5 6]

    # Coordinates
    coord = [1 -1.0 0.0 
        2 -4.0 0.0 
        3 -4.0 4.0 
        4 -0.7071 0.7071 
        5 0.0 4.0 
        6 0.0 1.0]

    return conn,coord
end