from simnibs import sim_struct

def prepare_emp(exp, tms=False):
    S = sim_struct.SESSION()
    if tms:
        sens_list = S.add_tmslist()
        add_func = sens_list.add_position
    else:
        sens_list = S.add_tdcslist()
        add_func = sens_list.add_electrode
    if exp == "P1":
        sens_list.currents = [1e-3, -1e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'P8'
        elec1.pos_ydir = 'TP8'
        elec1.shape = 'rect'
        elec1.dimensions = [70, 50]
        elec1.thickness = 5
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF7'
        elec2.pos_ydir = 'Fpz'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P2":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'P3'
        elec1.pos_ydir = 'PO7'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 70]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'Rch'
        elec2.pos_ydir = 'yRch'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P2-5050":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'P3'
        elec1.pos_ydir = 'PO7'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 50]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'Rch'
        elec2.pos_ydir = 'yRch'
        elec2.shape = 'rect'
        elec2.dimensions = [50, 50]
        elec2.thickness = 5
    elif exp == "P3":
        sens_list.currents = [1e-3, -1e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'CP5'
        elec1.pos_ydir = 'C5'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 70]
        elec1.thickness = 5
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF4'
        elec2.pos_ydir = 'Fp2'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P4":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'FC5'
        elec1.pos_ydir = 'F5'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 70]
        elec1.thickness = 5
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF4'
        elec2.pos_ydir = 'Fp2'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P5":
        sens_list.currents = [1e-3, -1e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'C3'
        elec1.pos_ydir = 'FC3'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 70]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF4'
        elec2.pos_ydir = 'Fp2'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P5-2ma":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'C3'
        elec1.pos_ydir = 'FC3'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 70]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF4'
        elec2.pos_ydir = 'Fp2'
        elec2.shape = 'rect'
        elec2.dimensions = [70, 50]
        elec2.thickness = 5
    elif exp == "P5-2ma-5050":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'C3'
        elec1.pos_ydir = 'FC3'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 50]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'AF4'
        elec2.pos_ydir = 'Fp2'
        elec2.shape = 'rect'
        elec2.dimensions = [50, 50]
        elec2.thickness = 5
    elif exp == "P6":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'I2'
        elec1.pos_ydir = 'P4'
        elec1.shape = 'rect'
        elec1.dimensions = [50, 50]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'Rbu'
        elec2.pos_ydir = 'yRbu'
        elec2.shape = 'rect'
        elec2.dimensions = [50, 50]
        elec2.thickness = 5
    elif exp == "P6-3030":
        sens_list.currents = [2e-3, -2e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'I2'
        elec1.pos_ydir = 'P4'
        elec1.shape = 'rect'
        elec1.dimensions = [30, 30]
        elec1.thickness = 5
        # second electrode
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'Rbu'
        elec2.pos_ydir = 'yRbu'
        elec2.shape = 'rect'
        elec2.dimensions = [30, 30]
        elec2.thickness = 5
    elif exp == "P7":
        sens_list.currents = [1e-3, -1e-3]
        # first electrode
        elec1 = add_func()
        elec1.channelnr = 1
        elec1.centre = 'F4'
        elec1.shape = 'ellipse'
        elec1.dimensions = [25, 25]
        elec1.thickness = 5
        elec2 = add_func()
        elec2.channelnr = 2
        elec2.centre = 'F4'
        elec2.shape = 'ellipse'
        elec2.dimensions = [115, 115]
        elec2.thickness = 5
        hole = elec2.add_hole()
        hole.centre = "F4"
        hole.shape = "ellipse"
        hole.dimensions = [92, 92]
    elif exp == "P8":
        centres = ["F3", "F1", "AF3", "F5", "FC3"]
        sens_list.currents = [2e-3, -.5e-3, -.5e-3, -.5e-3, -.5e-3]
        for idx, centre in zip(range(1,6), centres):
            elec = add_func()
            elec.channelnr = idx
            elec.shape = "ellipse"
            elec.dimensions = [10, 10]
            elec.thickness = 5
            elec.centre = centre

    return S
