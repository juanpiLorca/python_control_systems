import json 
import csv 
from cell import Cell

class Line(): 
    """
    Runs sequentially simulations in each cell (non-parallel) ==> it might be a problem 
    Args: 
        >>> num_cells: amount of flotation cells to be instantiate 
        >>> num_ore_classes: number of mineralogical classes
        >>> num_bubble_classes: number of bubble classes
        >>> params: Cell's parameters 
    """

    def __init__(self, num_cells: int, num_ore_classes: int, num_bubble_classes: int, 
                 params: dict): 
        self.num_cells = num_cells
        self.num_ores = num_ore_classes
        self.num_classes = num_bubble_classes

        # Cell's parameters: 
        self.params = params

        for i in range(0, self.num_cells): 
            if i == 0: 
                Qfeed = self.params["Qfeed"]
            else: 
                Qfeed = 0
            vars(self)[f"Cell_{i}"] = Cell(cell_id=i, i=self.num_ores, 
                                           num_classes=self.num_classes, 
                                           params=self.params, Qfeed=Qfeed) 

    def update(self): 
        Cell_0 = vars(self)["Cell_0"]
        x0 = Cell_0.simulate()

        Qtail_i = Cell_0.Qtail
        for i in range(1, self.num_cells): 
            Cell_i = vars(self)[f"Cell_{i}"]
            Cell_i.Qfeed = Qtail_i
            xi = Cell_i.simulate()
            # Pull new Qtail value and assign it to the next Qfeed: 
            Qtail_i = Cell_i.Qtail



# Simulation: ----------------------------------------------
NUM_CELLS = 3
fieldnames = [f"y{i}_value" for i in range(NUM_CELLS)]
fieldnames = ["time"] + fieldnames

with open("data_line.csv", "w") as csv_file:
    csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    csv_writer.writeheader() 

def main(num_cells=NUM_CELLS): 
    num_ore_classes = 2
    num_bubble_classes = 5
    with open("python-version/params.json") as fjson: 
        params = json.load(fjson)

    line = Line(num_cells=num_cells, num_ore_classes=num_ore_classes, 
                num_bubble_classes=num_bubble_classes, params=params)
    
    info = {f"y{i}_value": 0 for i in range(NUM_CELLS)}
    info["time"] = 0 

    while True: 
        line.update()
        t = vars(line)['Cell_0'].t
        info["time"] = t
        for i in range(NUM_CELLS): 
            info[f"y{i}_value"] = vars(line)[f"Cell_{i}"].x[0]

        with open("data_line.csv", "a") as csv_file: 
            csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames) 
            csv_writer.writerow(info)
        
        print(f"Opening valve 1: {round(vars(line)['Cell_0'].u[0], 5)} \
              Air flux 1: {round(vars(line)['Cell_0'].u[1], 5)}")
        print(f"Opening valve 2: {round(vars(line)['Cell_1'].u[0], 5)} \
              Air flux 2: {round(vars(line)['Cell_1'].u[1], 5)}")
        print(f"Opening valve 3: {round(vars(line)['Cell_2'].u[0], 5)} \
              Air flux 3: {round(vars(line)['Cell_2'].u[1], 5)}")
        #time.sleep(0.5)


if __name__ == "__main__": 
    main()