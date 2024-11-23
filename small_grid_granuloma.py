
from cc3d import CompuCellSetup
        

from small_grid_granulomaSteppables import small_grid_granulomaSteppable
CompuCellSetup.register_steppable(steppable=small_grid_granulomaSteppable(frequency=1))

from small_grid_granulomaSteppables import mac_recruitment_Steppable
CompuCellSetup.register_steppable(steppable=mac_recruitment_Steppable(frequency=3))

from small_grid_granulomaSteppables import InfectionSteppable
CompuCellSetup.register_steppable(steppable=InfectionSteppable(frequency=1))

from small_grid_granulomaSteppables import CellTypeConverterSteppable
CompuCellSetup.register_steppable(steppable=CellTypeConverterSteppable(frequency=1))

from small_grid_granulomaSteppables import inf_mac_killingSteppable
CompuCellSetup.register_steppable(steppable=inf_mac_killingSteppable(frequency=1))

from small_grid_granulomaSteppables import SurfaceSteppable
CompuCellSetup.register_steppable(steppable=SurfaceSteppable(frequency=1))

from small_grid_granulomaSteppables import Hypoxia_deathSteppable
CompuCellSetup.register_steppable(steppable=Hypoxia_deathSteppable(frequency=10))

from small_grid_granulomaSteppables import OutputFileSteppable
CompuCellSetup.register_steppable(steppable=OutputFileSteppable(frequency=50))

CompuCellSetup.run()
